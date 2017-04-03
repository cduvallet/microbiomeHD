# ${1} is the name of the target tar file, in format:
# data/raw_otu_tables/nash_chan_results.tar.gz
## Usage: copy_tar_folders.sh data/raw_otu_tables/target_folder.tar.gz

## This code will be replaced by a call to download from wherever
orig=${1#data/raw_otu_tables/}
orig=${orig%.tar.gz}
origdir=/Users/claire/github/disease_dataset_project/data/aws_processing_results

target=$1

tar -C $origdir -czvkf $target $orig

## This code will stay
targetdir=${target%/*}
tar -C $targetdir -xvkf $target
