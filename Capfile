
require 'catpaws'
#generic settings
set :aws_access_key,  ENV['AMAZON_ACCESS_KEY']
set :aws_secret_access_key , ENV['AMAZON_SECRET_ACCESS_KEY']
set :ec2_url, ENV['EC2_URL']
set :ssh_options, { :user => "ubuntu", :keys=>[ENV['EC2_KEYFILE']]}
set :key, ENV['EC2_KEY']
set :key_file, ENV['EC2_KEYFILE']
set :ami, '' 
set :instance_type, ''
set :working_dir, '/mnt/work'
set :nhosts, 1
set :group_name, ''
set :snap_id, `cat SNAPID`.chomp 
set :vol_id, `cat VOLUMEID`.chomp 
set :ebs_size, 10 
set :availability_zone, 'eu-west-1a'
set :dev, '/dev/sdf'
set :mount_point, '/mnt/data'


# Allow Capfile.local to override these settings
begin
 load("Capfile.local")
rescue Exception
end

#make a new EBS volume from this snap 
#cap EBS:create

#and mount your EBS
#cap EBS:attach
#cap EBS:mount_xfs

#if you want to keep the results

#cap EBS:snapshot


#and then shut everything down:

# cap EBS:unmount
# cap EBS:detach
# cap EBS:delete - unless you're planning to use it again.
# cap EC2:stop

    