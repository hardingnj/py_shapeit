__author__ = 'Nicholas Harding'

import re
from os.path import join, isfile, isdir
from os import mkdir, getenv
import sh
import uuid
import tempfile
from utils import parse_command, create_sh_script, md5_for_file
import subprocess
import yaml
import sys


class ShapeIt():

    @staticmethod
    def _get_version(executable):

        p = subprocess.Popen(args=[executable, '--version'],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        out, _ = p.communicate()

        p = re.compile("Version : (.+)\n")  # parentheses for capture groups
        m = p.search(out)
        if m:
            return m.group(1)
        else:
            print(out)
            raise Exception('Version not parsed.')

    def __init__(self, outdir, executable=getenv('SHAPEIT', 'shapeit'),
                 ligate_exec=getenv('LIGATEHAPLOTYPES', 'ligateHAPLOTYPES'),
                 duohmm_exec=getenv('DUOHMM', 'duohmm'),
                 version=None):

        self.executable = executable
        self.ligate_bin = ligate_exec
        self.duohmm_bin = duohmm_exec
        if version is None:
            version = self._get_version(self.executable)
        self.version = version
        self.name = 'ShapeIt'
        self.checksum_file = None
        self.contig = None

        self.run_id = self.name + '_' + str(uuid.uuid4().get_hex().upper()[0:8])

        self.outdir = join(outdir, self.name, self.version, self.run_id)
        self.basedir = outdir
        self.si_job_list = list()
        self.ligate_script = None
        self.duohmm_script = None
        self.h5_script = None
        self.si_script = None
        self.dirs = {d: join(self.outdir, d) for d in ('log', 'script')}
        [sh.mkdir(d, '-p') for d in self.dirs.values() if not isdir(d)]

        self.haplotypes_f = join(self.outdir, self.run_id + '.haps.gz')
        self.phased_f = join(self.outdir, self.run_id + '.samples.gz')

        self.h5out = join(self.outdir, self.run_id + '.h5')

        self.duohmm_haps = None
        self.duohmm_sample = None
        self.duohmm_geno_e = None
        self.duohmm_recomb = None
        self.provided_sample_file = None

        self.param_f = join(self.outdir, self.run_id + '.yaml')
        self.settings = {'run_id': self.run_id, 'version': self.version,
                         'executable': self.executable, 'outdir': self.outdir,
                         'basedir': outdir, 'haps': self.haplotypes_f,
                         'samples': self.phased_f}

    def setup_region_jobs(self, parameters, regions=None, contig=None,
                          vcf_file=None, pirs=None, graphs=False,
                          duohmm=False, sample_file=None):

        """
        Setup the directory structure and scripts for a shapeIt job that
          automatically splits the job into regions and optionally uses PIRs
          and duoHMM.
        :param parameters: List of parameters to be passed to shapeIt. See
         SHAPEIT docs for full details.
        :param regions: A list of tuples of type (start, end) that represents
         the regions the job is split into. Should have at least 5kb overlap,
         so ligateHAPLOTYPES can work
        :param contig: An ID for the contig being phased
        :param vcf_file:
          VCF format 4.0 file. In gzipped or plaintext. IMPORTANT: SHAPEIT does
          not support bzgipped format VCF files.
        :param pirs:
          A string representing PIR files that contain the phase informative
          reads for each region in the region list. The string should be
          formatted with {0}{1} where 0/1 are replaced with start/end of the
          region.
        :param graphs:
          boolean, whether to produce haplotype graphs for each region.
        :param duohmm:
          boolean, whether to use the duoHMM tool to correct haplotypes using
          pedigree information.
        :param sample_file:
          if duohmm, the sample file containing the pedigree information.
        :return:
        """

        region_dir = join(self.outdir, 'region_outputs')
        mkdir(region_dir)

        parameters = [str(x) for x in parameters] + ['--input-vcf', vcf_file]
        self.checksum_file = vcf_file
        self.provided_sample_file = sample_file
        self.contig = contig

        hap_files = list()
        for i, region in enumerate(regions):
            start, stop = [str(x) for x in region]
            region_stem = str(i) + '_' + self.run_id

            pir_region = None
            if pirs is not None:
                pir_region = pirs.format(start, stop)

            script_name, haps, _ = self._setup_shapeit(outdir=region_dir,
                                                       outstem=region_stem,
                                                       si_parameters=parameters,
                                                       start=start, stop=stop,
                                                       pirfile=pir_region,
                                                       graphs=graphs)

            # name, script, mem, dependency
            hap_files.append(haps)
            self.si_job_list.append(script_name)

        # set up ligateHaplotypes
        self._setup_ligate_haplotypes(hap_files, vcf_file)

        # set up duohmm if specified
        if duohmm:
            self._setup_duohmm()

        # run shapeIt to hdf5 script
        self._setup_h5_conversion()

        # finally add parameters to settings doc
        self.settings['params'] = parse_command(parameters)

    def run_region_qsub(self, si_args, li_args, dm_args=None):
        """
        Run the previously created jobs via qsub. This is designed for use with
        SGE, perhaps a later version will have more general dependency support.
        :param si_args: Qsub arguments passed to shapeIt jobs.
        :param li_args: Qsub arguments passed to the ligateHAPLOTYPE job.
        :param dm_args: Qsub arguments passed to the duohmm job.
        :return: None
        """
        qsub_parameters = ['-S', '/bin/bash',
                           '-j', 'y',
                           '-o', self.dirs['log']]
        hold_list = [self.run_id]

        # get checksum
        self.settings['checksum'] = md5_for_file(self.checksum_file)
        yaml.dump(self.settings,
                  stream=open(self.param_f, 'w'))

        for job in self.si_job_list:
            print sh.qsub('-N', self.run_id, qsub_parameters, si_args, job)

        print sh.qsub('-hold_jid', self.run_id, '-N', 'ligate' + self.run_id,
                      qsub_parameters + li_args, self.ligate_script)
        hold_list.append('ligate' + self.run_id)

        if self.duohmm_script is not None:
            print sh.qsub('-hold_jid', ",".join(hold_list),
                          '-N', 'duohmm' + self.run_id,
                          qsub_parameters + dm_args, self.duohmm_script)
            hold_list.append('duohmm' + self.run_id)

        print sh.qsub('-hold_jid', ",".join(hold_list),
                      '-N', 'si2hdf5' + self.run_id, '-l', 'h_vmem=4G',
                      qsub_parameters, self.h5_script)

    def setup_single_job(self, parameters, contig, vcf_file,
                         start=None, stop=None, pirs=None,
                         graphs=False, duohmm=False, sample_file=None):

        # basically, set up a shapeIT job running the whole file as one.
        parameters = [str(x) for x in parameters] + ['--input-vcf', vcf_file]
        self.checksum_file = vcf_file
        self.provided_sample_file = sample_file
        self.contig = contig

        script, haps, sample = self._setup_shapeit(self.outdir, self.run_id,
                                                   parameters, start=start,
                                                   stop=stop, pirfile=pirs,
                                                   graphs=graphs)
        self.si_script = script
        assert haps == self.haplotypes_f
        assert sample == self.phased_f

        if duohmm:
            self._setup_duohmm()

        self._setup_h5_conversion()

        self.settings['params'] = parse_command(parameters)

    def run_single_qsub(self, si_args, dm_args=None):

        hold_list = list()
        qsub_parameters = ['-S', '/bin/bash',
                           '-j', 'y',
                           '-o', self.dirs['log']]
        # get checksum
        self.settings['checksum'] = md5_for_file(self.checksum_file)
        yaml.dump(self.settings,
                  stream=open(self.param_f, 'w'))

        print sh.qsub('-N', self.run_id, qsub_parameters,
                      si_args, self.si_script)
        hold_list.append(self.run_id)

        if self.duohmm_script is not None:
            print sh.qsub('-hold_jid', self.run_id,
                          '-N', 'duohmm' + self.run_id,
                          qsub_parameters + dm_args, self.duohmm_script)
            hold_list.append('duohmm'+self.run_id)

        print sh.qsub('-hold_jid', ",".join(hold_list),
                      '-N', 'si2hdf5' + self.run_id, '-l', 'h_vmem=4G',
                      qsub_parameters, self.h5_script)

    # INTERNAL FUNCTIONS
    def _setup_duohmm(self):

        gunzip = "gunzip -c {0} > {1}"
        tmpduo = tempfile.NamedTemporaryFile(delete=False)
        self.duohmm_script = join(self.dirs['script'], 'duohmm.sh')

        duohmm_root = join(self.outdir, self.run_id + '_duohmm')
        self.duohmm_haps = duohmm_root + '.haps'
        self.duohmm_sample = duohmm_root + '.samples'
        self.duohmm_geno_e = duohmm_root + '.GE.txt'
        self.duohmm_recomb = duohmm_root + '.RC.txt'

        cmd_duohmm = [self.duohmm_bin,
                      '-H', tmpduo.name,
                      '-O', duohmm_root,
                      '-G', self.duohmm_geno_e,
                      '-R', self.duohmm_recomb]

        command_list = [gunzip.format(self.haplotypes_f, tmpduo.name + '.haps'),
                        "cp {0} {1}".format(self.provided_sample_file,
                                            tmpduo.name + '.sample'),
                        'cd ' + self.outdir,
                        " ".join(cmd_duohmm),
                        " ".join(["gzip", self.duohmm_haps, self.duohmm_sample,
                                  self.duohmm_geno_e, self.duohmm_recomb])]

        create_sh_script(self.duohmm_script, command_list,
                         outfile=self.duohmm_haps + '.gz')

    def _setup_h5_conversion(self):

        self.h5_script = join(self.dirs['script'], 'convert_h5.sh')
        h5_root = "{0} {1}/bin/shapeIt2hdf5.py {2} {3} --chr {4} --out {5}"

        if self.duohmm_script is not None:
            haps = self.duohmm_haps + '.gz'
            samp = self.duohmm_sample + '.gz'
        else:
            haps = self.haplotypes_f
            samp = self.phased_f

        h5_cmd = h5_root.format(sys.executable, sys.exec_prefix, haps, samp,
                                self.contig, self.h5out)

        create_sh_script(self.h5_script, [h5_cmd], self.h5out)

    def _setup_ligate_haplotypes(self, haplotype_files, vcf):

        tmp = tempfile.NamedTemporaryFile(delete=False)

        self.ligate_script = join(self.dirs['script'], 'ligatehaps.sh')

        cmd_ligate = [self.ligate_bin, '--vcf', tmp.name, '--chunks',
                      " ".join(haplotype_files), '--output', self.haplotypes_f,
                      self.phased_f]

        # ligateHAPLOTYPES only accepts uncompressed VCF files
        gunzip = "gunzip -c {0} > {1}"
        create_sh_script(self.ligate_script,
                         [gunzip.format(vcf, tmp.name), 'cd ' + self.outdir,
                          " ".join(cmd_ligate), 'rm ' + tmp.name],
                         outfile=self.haplotypes_f)

    def _setup_shapeit(self, outdir, outstem, si_parameters, start=None,
                       stop=None, pirfile=None, graphs=False):

        # this is pretty much the core function
        haps = join(outdir, outstem + '.haps.gz')
        samples = join(outdir, outstem + '.samples.gz')

        cmd_shape_it = [self.executable]

        # below code also applies to single
        if pirfile is not None:
            if graphs:
                raise ValueError("PIRs not compatible with graphs")
            assert isfile(pirfile)
            cmd_shape_it += ['-assemble', '--input-pir', pirfile]

        if (start is not None) or (stop is not None):
            assert (stop is not None) and (start is not None), \
                "If start/stop provided other must be as well."
            cmd_shape_it += ['--input-from', start, '--input-to', stop]

        cmd_shape_it += si_parameters

        if graphs:
            graph_fn = outstem + '.graph'
            cmd_shape_it += ['--output-graph', graph_fn]

            graph_convert = "{shapeit} -convert --input-graph {input} " + \
                "--output-max {phased} {sample} --aligned"

            convert_cmd = graph_convert.format(shapeit=self.executable,
                                               input=graph_fn,
                                               phased=haps,
                                               sample=samples)
        else:
            cmd_shape_it += ['--output-max', haps, samples]
            convert_cmd = ''

        script_name = join(self.dirs['script'], outstem + '_shapeIt.sh')

        exit_check = "\n".join(["{",
                                "if [ -f %s.ok ]; then",
                                "\techo \"Already completed ok!\"",
                                "\texit 0",
                                "fi",
                                "}"])

        create_sh_script(
            filename=script_name,
            commands=['cd ' + outdir, exit_check % haps,
                      " ".join(cmd_shape_it), convert_cmd],
            outfile=haps)

        return script_name, haps, samples