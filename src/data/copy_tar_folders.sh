# ${1} is the name of the target tar file, in format:
# data/raw_otu_tables/nash_chan_results.tar.gz

orig=${1#data/raw_otu_tables/}
orig=${orig%.tar.gz}
origdir=/Users/claire/github/disease_dataset_project/data/aws_processing_results

target=$1

tar -C $origdir -czvkf $target $orig
