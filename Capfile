###
# Align the data with Bowtie and dump the resulting BAM files
# on S3
# Run Macs peak finding on them.

require 'catpaws'

#generic settings
set :aws_access_key,  ENV['AMAZON_ACCESS_KEY']
set :aws_secret_access_key , ENV['AMAZON_SECRET_ACCESS_KEY']
set :ec2_url, ENV['EC2_URL']
set :ssh_options, { :user => "ubuntu", :keys=>[ENV['EC2_KEYFILE']]}
set :key, ENV['EC2_KEY']
set :key_file, ENV['EC2_KEYFILE']
set :ami, `curl http://mng.iop.kcl.ac.uk/cass_data/buckley_ami/AMIID`.chomp
set :instance_type, 'm2.4xlarge'
set :s3cfg, ENV['S3CFG'] #location of ubuntu s3cfg file
set :working_dir, '/mnt/work'
set :availability_zone, 'eu-west-1a'

#note to self
#ami-52794c26 32bit Lucid
#ami-505c6924 64bit Maverick
#ami-20794c54 64bit Lucid

set :nhosts, 1
set :group_name, 'ctx12_epigenetics'

set :snap_id, `cat SNAPID`.chomp #ec2 eu-west-1 
set :vol_id, `cat VOLUMEID`.chomp #empty until you've created a new volume
set :ebs_size, 100  #Needs to be the size of the snap plus enough space for alignments
set :ebs_zone, 'eu-west-1a'  #is where the ubuntu ami is
#set :dev, '/dev/sdh'
#######comment out the sdh line after creating and attaching the EC2 and EBS and before formatting and mounting###########
set :dev, '/dev/xvdh'
set :mount_point, '/mnt/data'

###ssh to server
#

###sft to server
#sftp -o"IdentityFile=/home/spkaabi/ec2/spkaabi.pem" ubuntu@ec2-46-51-154-242.eu-west-1.compute.amazonaws.com

#start EC2 instances
#cap EC2:start

#make a new EBS volume from this snap 
#cap EBS:create

#and attach your EBS
#cap EBS:attach

#format a new EBS
#cap EBS:format_xfs

#and mount your EBS
#cap EBS:mount_xfs

#### Data uploading.

desc "Upload data files"
task :upload_data, :roles => group_name do
    run "rsync -e 'ssh -i /home/spkaabi/ec2/spkaabi.pem' -vzP data/* ubuntu@ec2-46-51-154-242.eu-west-1.compute.amazonaws.com:#{mount_point}" 
end
before 'upload_data', 'EC2:start'

desc "unzip data"
task :unzip_data, :roles => group_name do
    files = capture("ls #{mount_point}/*.gz")
  files = files.split("\n")
  files.each {|f| 
      run "cd #{mount_point}/ && gunzip #{f}"
  }
end
before 'unzip_data', 'EC2:start' 

#### Quality Control

# convert the export.txt file to fastq for alignment
task :make_fastq, :roles => group_name do 
  upload("/space/cassj/chipseq_pipeline/export2fastq.pl", "#{working_dir}/export2fastq.pl")
  run "chmod +x #{working_dir}/export2fastq.pl"
  run "sudo mv #{working_dir}/export2fastq.pl /usr/local/bin/"

  files = capture("ls #{mount_point}/*export.txt").split("\n")
  files = files.map {|f| f.chomp}

  files.each{|infile| 
    outfile = infile.sub('.txt', '.fastq')
    run "export2fastq.pl #{infile} > #{outfile}"
  } 

end
before 'make_fastq', 'EC2:start' 

# Run fastQC on the data files
desc "run fastqc"
task :fastqc, :roles => group_name do
  files = capture("ls #{mount_point}/*.fq").split("\n")
  files = files.map {|f| f.chomp}
   
  files.each{|infile| 
    run "fastqc --outdir #{mount_point} #{infile}"
  } 

end
before 'fastqc', 'EC2:start'

# Pull the results back to the mng.iop.kcl.ac.uk server
desc "download fastqc files"
task :get_fastqc, :roles => group_name do
  `rm -Rf results/fastqc` #remove previous results
  `mkdir -p results/fastqc`
  files = capture "ls #{mount_point}/*fastqc.zip"
  files = files.split("\n")
  files.each{|f|
    outfile = f.sub(/.*\//,'')
    download( "#{f}", "results/fastqc/#{outfile}")
    `cd results/fastqc && unzip #{outfile} && rm #{outfile}`
  }
end
before "get_fastqc", 'EC2:start'

#### Alignment


#get the current human genome from bowtie prebuilt indexes
task :fetch_human_genome, :roles => group_name do
  run "mkdir -p #{working_dir}/indexes"
  run "cd #{working_dir}/indexes && curl ftp://ftp.cbcb.umd.edu/pub/data/bowtie_indexes/hg19.ebwt.zip > hg19.ebwt.zip"
  run "rm -Rf #{working_dir}/indexes/chr*"
  run "cd  #{working_dir}/indexes && unzip -o hg19.ebwt.zip"
#?  run "export BOWTIE_INDEXES='#{working_dir}/indexes'"
end
before "fetch_human_genome","EC2:start"

# s3 setup
desc "install s3 client"
task :install_s3, :roles => group_name do
  sudo 'apt-get update'
  sudo 'apt-get install -y s3cmd'
end
before 'install_s3', 'EC2:start'

desc "upload the s3 config file"
task :s3_config, :roles => group_name do
  upload("/home/kkvi1130/.s3cfg", '/home/ubuntu/.s3cfg')
end 
before 's3_config', 'EC2:start'


#get the current mouse genome (which I already have on S3).
task :fetch_mouse_genome, :roles => group_name do
  run "mkdir -p #{working_dir}/indexes"
  run "cd #{working_dir}/indexes && curl https://s3-eu-west-1.amazonaws.com/genome-mm9-bowtie/mm9.ebwt.zip > mm9.ebwt.zip"
  run "rm -Rf #{working_dir}/indexes/chr*"
  run "cd  #{working_dir}/indexes && unzip -o mm9.ebwt.zip"
end 
before "fetch_mouse_genome","EC2:start"

# run bowtie on the fastq file
# This is recent illumina data, quals should be post v1.3
task :run_bowtie, :roles => group_name do

  files = capture("ls #{mount_point}/*.fq").split("\n")
  files = files.map {|f| f.chomp}

  files.each{|infile|
    outfile = infile.sub('.fq', '.sam')
    run ("export BOWTIE_INDEXES='#{working_dir}/indexes' && bowtie  --sam --best -k1 -l15 -n1 -m3 -p20 --solexa1.3-quals --chunkmbs 256  -q mm9 --quiet  #{infile}  > #{outfile} ")
  } 

end
before "run_bowtie", "EC2:start"

# Make binary BAM files from SAM
desc "make bam from sam"
task :to_human_bam, :roles => group_name do
  run "curl 'http://github.com/cassj/my_bioinfo_scripts/raw/master/genomes/hg19_lengths' > #{working_dir}/hg19_lengths"
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/\.sam$/)}
  files.each{|f| 
    f_out = f.sub('.sam', '.bam')
    run "samtools view -bt #{working_dir}/hg19_lengths -o #{mount_point}/#{f_out} #{mount_point}/#{f}"
  }
end
before "to_human_bam", "EC2:start"

desc "make bam from sam"
task :to_mouse_bam, :roles => group_name do
  run "curl 'http://github.com/cassj/my_bioinfo_scripts/raw/master/genomes/mm9_lengths' > #{working_dir}/mm9_lengths"
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/\.sam$/)}
  files.each{|f| 
    f_out = f.sub('.sam', '.bam')
    run "samtools view -bt #{working_dir}/mm9_lengths -o #{mount_point}/#{f_out} #{mount_point}/#{f}"
  }
end
before "to_mouse_bam", "EC2:start"

# Sort the BAM files
desc "sort bam"
task :sort_bam, :roles => group_name do
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/\.bam/)}
  files.each{|f| 
    f_out = f.sub('.bam', '_sorted')
    run "cd #{mount_point} && samtools sort #{f}  #{f_out}"
  }
end
before "sort_bam", "EC2:start"


# Remove PCR Duplicate Reads
desc "remove duplicates"
task :rmdups, :roles => group_name do
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted\.bam/)}
  files.each{|f| 
    f_out = f.sub('sorted', 'sorted_nodups')
    run "cd #{mount_point} && samtools rmdup -s #{f}  #{f_out}"
  }
end
before "rmdups", "EC2:start"



# Index the BAM files
desc "index bam files"
task :index, :roles => group_name do
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted_nodups\.bam/)}
  files.each{|f| 
    f_out = f.sub('.bam', '.bai')
    run "cd #{mount_point} && samtools index  #{f} #{f_out}"
  }
end
before "index", "EC2:start"

# Create a summary of the files
desc "create a summary of the bam files"
task :flagstat, :roles => group_name do
 files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted_nodups\.bam/)}
  files.each{|f|
    f_out = f.sub('.bam', '.summary')
    run "cd #{mount_point} && samtools flagstat #{f} > #{f_out}"
  }

end
before "flagstat", "EC2:start"


# Pull the BAM files back to the mng.iop.kcl.ac.uk server
desc "download bam files"
task :get_bam, :roles => group_name do
  `rm -Rf results/alignment/bowtie` #remove previous results
  `mkdir -p results/alignment/bowtie`
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted_nodups/)}
  files.each{|f|
    download( "#{mount_point}/#{f}", "results/alignment/bowtie/#{f}")
  }
end
before "get_bam", "EC2:start"

##so don't need to do one by one - after make bam files
#before "get_bam","EC2:start","sort_bam","rmdups","index","flagstat","EBS:snapshot"


### Macs ?

#macs_url ="http://liulab.dfci.harvard.edu/MACS/src/MACS-1.4.0beta.tar.gz"
#macs_version = "MACS-1.4.0beta"
#
#task :install_macs, :roles => group_name do
#  sudo "apt-get install -y python"
#  run "cd #{working_dir} && wget --http-user macs --http-passwd chipseq #{macs_url}"
#  run "cd #{working_dir} && tar -xvzf #{macs_version}.tar.gz"
#  run "cd #{working_dir}/#{macs_version} && sudo python setup.py install"
#  sudo "ln -s /usr/local/bin/macs* /usr/local/bin/macs"
#end
#before "install_macs", 'EC2:start'
#
#task :install_peaksplitter, :roles => group_name do
#  url ='http://www.ebi.ac.uk/bertone/software/PeakSplitter_Cpp_1.0.tar.gz'
#  filename = 'PeakSplitter_Cpp_1.0.tar.gz'
#  bin = 'PeakSplitter_Cpp/PeakSplitter_Linux64/PeakSplitter'
#  run "cd #{working_dir} && curl #{url} > #{filename}"
#  run "cd #{working_dir} && tar -xvzf #{filename}"
#  run "sudo cp #{working_dir}/#{bin} /usr/local/bin/PeakSplitter"
#end 
#before 'install_peaksplitter', 'EC2:start'

desc "run macs"
task :run_macs, :roles => group_name do
  ChIP = "#{mount_point}/REST_Astro.clean_sorted_nodups.bam"
  Input = "#{mount_point}/input_Astro.clean_sorted_nodups.bam"
  
  genome = 'mm'
  bw = [300]
  

      dir = "#{mount_point}/macs"
      run "rm -Rf #{dir}"
      run "mkdir #{dir}"

      macs_cmd =  "macs --treatment #{ChIP} --control #{Input} --name REST_NS5dAstro --format BAM --gsize #{genome}"
      run "cd #{dir} && #{macs_cmd}"
  
end
before 'run_macs', 'EC2:start'


desc "bamToBed"
task :bamToBed, :roles => group_name do
   files = capture("ls #{mount_point}/*sorted_nodups.bam").split("\n")
   files = files.map {|f| f.chomp}
   files.each{|infile|
     f_out = infile.sub('.bam', '.bed')
     run "bamToBed -i #{infile} > #{f_out}"
   }

end
before 'bamToBed', 'EC2:start'

####run SICER
desc "run SICER"
task :run_SICER, :roles => group_name do

   ctx12_input    = 'CTX_input.clean_sorted_nodups.bed'
   bmp4_input     = 'BMP_input.clean_sorted_nodups.bed'
   fbs_input      = 'FBS_input.clean_sorted_nodups.bed'

   ctx12_h3k4me3  = 'CTX_H3K4.clean_sorted_nodups.bed'
   bmp4_h3k4me3   = 'BMP_H3K4.clean_sorted_nodups.bed'
   fbs_h3k4me3    = 'FBS_H3K4.clean_sorted_nodups.bed'

   ctx12_h3k27me3 = 'CTX_H3K27.clean_sorted_nodups.bed'
   bmp4_h3k27me3  = 'BMP_H3K27.clean_sorted_nodups.bed'
   fbs_h3k27me3   = 'FBS_H3K4.clean_sorted_nodups.bed'

   species = 'mm9'
   thresh = 1
   window_size = 200
   fragment_size = 300
   effective_genome_fraction = '0.75' 
   gap_size_k4 = 200
   gap_size_k27 = 600
   FDR = '0.001'

# /usr/local/bin/SICER [InputDir] [bed file] [control file] [OutputDir] [Species] [redundancy threshold] [window size (bp)] [fragment size] [effective genome fraction] [gap size (bp)] [FDR]

   run "mkdir -p #{mount_point}/SICER"


   run "mkdir -p #{mount_point}/SICER/ctx12_h3k4me3" 
   run "cd #{mount_point}/SICER/ctx12_h3k4me3 && SICER #{mount_point} #{ctx12_h3k4me3} #{ctx12_input} #{mount_point}/SICER/ctx12_h3k4me3 #{species} #{thresh} #{window_size} #{fragment_size} #{effective_genome_fraction} #{gap_size_k4} #{FDR}"

   run "mkdir -p #{mount_point}/SICER/ctx12_h3k27me3"
   run "cd #{mount_point}/SICER/ctx12_h3k27me3 && SICER #{mount_point} #{ctx12_h3k27me3} #{ctx12_input} #{mount_point}/SICER/ctx12_h3k27me3 #{species} #{thresh} #{window_size} #{fragment_size} #{effective_genome_fraction} #{gap_size_k27} #{FDR}"

   run "mkdir -p #{mount_point}/SICER/bmp4_h3k4me3"
   run "cd #{mount_point}/SICER/bmp4_h3k4me3 && SICER #{mount_point} #{bmp4_h3k4me3} #{bmp4_input} #{mount_point}/SICER/bmp4_h3k4me3 #{species} #{thresh} #{window_size} #{fragment_size} #{effective_genome_fraction} #{gap_size_k4} #{FDR}"

   run "mkdir -p #{mount_point}/SICER/bmp4_h3k27me3"
   run "cd #{mount_point}/SICER/bmp4_h3k27me3 && SICER #{mount_point} #{bmp4_h3k27me3} #{bmp4_input} #{mount_point}/SICER/bmp4_h3k27me3 #{species} #{thresh} #{window_size} #{fragment_size} #{effective_genome_fraction} #{gap_size_k27} #{FDR}"

   run "mkdir -p #{mount_point}/SICER/fbs_h3k4me3"
   run "cd #{mount_point}/SICER/fbs_h3k4me3 && SICER #{mount_point} #{fbs_h3k4me3} #{fbs_input} #{mount_point}/SICER/fbs_h3k4me3 #{species} #{thresh} #{window_size} #{fragment_size} #{effective_genome_fraction} #{gap_size_k4} #{FDR}"

   run "mkdir -p #{mount_point}/SICER/fbs_h3k27me3"
   run "cd #{mount_point}/SICER/fbs_h3k27me3 && SICER #{mount_point} #{fbs_h3k27me3} #{fbs_input} #{mount_point}/SICER/fbs_h3k27me3 #{species} #{thresh} #{window_size} #{fragment_size} #{effective_genome_fraction} #{gap_size_k27} #{FDR}"


end
before 'run_SICER', 'EC2:start'

#SICER /mnt/data CME141_GA3R71_export_sorted_nodups.bed CME140_s_5_export_sorted_nodups.bed /mnt/data/SICER/ctrl_h3k9ac mm9 1 200 300 0.75 600 0.001



#pack up the runs and downloads them to the server (without the wig files)
#task :pack_macs, :roles => group_name do
 # macs_dirs = capture "ls #{mount_point}"
  #macs_dirs = macs_dirs.split("\n").select {|f| f.match(/.*macs.*/)}
  #macs_dirs.each{|d|
   # run "cd #{mount_point} &&  tar --exclude *_wiggle* -cvzf #{d}.tgz #{d}"
  #}
  
#end
#before 'pack_macs','EC2:start' 

#task :get_macs, :roles => group_name do
 # macs_files = capture "ls #{mount_point}"
  #macs_files = macs_files.split("\n").select {|f| f.match(/.*macs.*\.tgz/)}
  #res_dir = 'results/alignment/bowtie/peakfinding/macs'
  #`rm -Rf #{res_dir}`
  #`mkdir -p #{res_dir}`
  #macs_files.each{|f| 
  #  download("#{mount_point}/#{f}", "#{res_dir}/#{f}") 
  #  `cd #{res_dir} && tar -xvzf #{f}`
  #}

#end
#before 'get_macs', 'EC2:start'


desc "convert peaks to IRanges"
task :peaks_to_iranges, :roles => group_name do
  run "cd #{working_dir} && rm -f peaksBed2IRanges_SICER.R"
  upload('scripts/peaksBed2IRanges_SICER.R' , "#{working_dir}/peaksBed2IRanges_SICER.R")
  run "cd #{working_dir} && chmod +x peaksBed2IRanges_SICER.R"

  run "cd #{mount_point}/SICER/ctx12_h3k4me3 && Rscript #{working_dir}/peaksBed2IRanges_SICER.R *islands-summary-FDR0.001"
  run "cd #{mount_point}/SICER/ctx12_h3k27me3 && Rscript #{working_dir}/peaksBed2IRanges_SICER.R *islands-summary-FDR0.001"
  run "cd #{mount_point}/SICER/bmp4_h3K4me3 && Rscript #{working_dir}/peaksBed2IRanges_SICER.R *islands-summary-FDR0.001"
  run "cd #{mount_point}/SICER/bmp4_h3k27me3 && Rscript #{working_dir}/peaksBed2IRanges_SICER.R *islands-summary-FDR0.001"
  run "cd #{mount_point}/SICER/fbs_h3k4me3 && Rscript #{working_dir}/peaksBed2IRanges_SICER.R *islands-summary-FDR0.001"
  run "cd #{mount_point}/SICER/fbs_h3k27me3 && Rscript #{working_dir}/peaksBed2IRanges_SICER.R *islands-summary-FDR0.001"

end
before 'peaks_to_iranges', 'EC2:start'

desc "annotate IRanges"
task :annotate_peaks, :roles => group_name do
  run "cd #{working_dir} && rm -f mm9RDtoGenes_toTSSfromSICER.R"
  upload('scripts/mm9RDtoGenes_toTSSfromSICER.R', "#{working_dir}/mm9RDtoGenes_toTSSfromSICER.R")
  run "cd #{working_dir} && chmod +x mm9RDtoGenes_toTSSfromSICER.R"
  
  run "cd #{mount_point}/SICER/ctx12_h3k4me3 && Rscript #{working_dir}/mm9RDtoGenes_toTSSfromSICER.R SICER.RangedData.RData"
  run "cd #{mount_point}/SICER/ctx12_h3k27me3 && Rscript #{working_dir}/mm9RDtoGenes_toTSSfromSICER.R SICER.RangedData.RData"
  run "cd #{mount_point}/SICER/bmp4_h3k4me3 && Rscript #{working_dir}/mm9RDtoGenes_toTSSfromSICER.R SICER.RangedData.RData"
  run "cd #{mount_point}/SICER/bmp4_h3k4me3 && Rscript #{working_dir}/mm9RDtoGenes_toTSSfromSICER.R SICER.RangedData.RData"
  run "cd #{mount_point}/SICER/fbs_h3k4me3 && Rscript #{working_dir}/mm9RDtoGenes_toTSSfromSICER.R SICER.RangedData.RData"
  run "cd #{mount_point}/SICER/fbs_h3k4me3 && Rscript #{working_dir}/mm9RDtoGenes_toTSSfromSICER.R SICER.RangedData.RData"


end
before 'annotate_peaks', 'EC2:start'

task :pack_peaks, :roles => group_name do
  macs_dirs = capture "ls #{mount_point}"
  macs_dirs = macs_dirs.split("\n").select {|f| f.match(/.*macs.*/)}
  macs_dirs.each{|d|
    unless d.match('tgz')
      run "cd #{mount_point} &&  tar --exclude *_wiggle* -cvzf #{d}.tgz #{d}"
    end
  }
  
end
before 'pack_peaks','EC2:start' 

task :get_peaks, :roles => group_name do
  macs_files = capture "ls #{mount_point}"
  macs_files = macs_files.split("\n").select {|f| f.match(/.*macs.*\.tgz/)}
  res_dir = 'results/alignment/bowtie/peakfinding/macs'
  `rm -Rf #{res_dir}`
  `mkdir -p #{res_dir}`
  macs_files.each{|f| 
    download("#{mount_point}/#{f}", "#{res_dir}/#{f}") 
    `cd #{res_dir} && tar -xvzf #{f}`
  }

end
before 'get_peaks', 'EC2:start'

# Create a summary of the files
desc "create a summary of the bam files"
task :flagstat, :roles => group_name do
 files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted_nodups\.bam/)}
  files.each{|f|
    f_out = f.sub('.bam', '.summary')
    run "cd #{mount_point} && samtools flagstat #{f} > #{f.out}"
  }

end
before "flagstat", "EC2:start"

task :peak_count, :roles => group_name do

  run "cd #{working_dir} && rm -f peak_read_counts.R"
  upload('scripts/peak_read_counts.R', "#{working_dir}/peak_read_counts.R")
  run "cd #{working_dir} && chmod +x peak_read_counts.R"

	rest_ip    = '#{mount_point}/CME143_GA3R71_export_sorted_nodups.bam'
	rest_peaks = '#{mount_point}/DeSeq_K9/rest_peaks_for_deseq.bed'
	ctrl_ip    = '#{mount_point}CME141_GA3R71_export_sorted_nodups.bam'
	ctrl_peaks = '#{mount_point}/DeSeq_K9/rest_peaks_for_deseq.bed'

  run "Rscript #{working_dir}/peak_read_counts.R #{mount_point}/DeSeq_K9 4 #{rest_ip} #{rest_peaks} #{ctrl_ip} #{ctrl_peaks}"

end
before 'peak_count', 'EC2:start'

task :compare_peaks, :roles => group_name do

### need to install edgeR first as sudo:
# in R source("http://bioconductor.org/biocLite.R")
#      biocLite("edgeR")

  run "cd #{working_dir} && rm -f compare_peaks_middlepeak.R"
  upload('scripts/compare_peaks.R', "#{working_dir}/compare_peaks_middlepeak.R")
  run "cd #{working_dir} && chmod +x compare_peaks_middlepeak.R"

  run  "Rscript #{working_dir}/compare_peaks_middlepeak.R #{mount_point}/DeSeq_K9/counts_middlepeak1000.RData #{mount_point}/DeSeq_K9/libsizes_middlepeak1000.RData #{mount_point}/DeSeq_K9"

end
before 'compare_peaks', 'EC2:start'





#if you want to keep the results
#cap EBS:snapshot


#and then shut everything down:

# cap EBS:unmount
# cap EBS:detach
# cap EBS:delete - unless you're planning to use it again.
# cap EC2:stop



