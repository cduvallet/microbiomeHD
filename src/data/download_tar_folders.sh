# ${1} is the name of the target tar file, in format:
# data/raw_otu_tables/nash_chan_results.tar.gz
## Usage: download_tar_folders.sh data/raw_otu_tables/target_folder.tar.gz

## This code will be replaced by a call to download from wherever
# If ${1} is data/raw_otu_tables/nash_chan_results.tar.gz, then orig is
# nash_chan_results.tar.gz
orig=${1#data/raw_otu_tables/}

# $1 is the full path, e.g. data/raw_otu_tables/nash_chan_results.tar.gz
target=$1

## Download nash_chan_results.tar.gz from Zenodo to data/raw_otu_tables/nash_chan_results.tar.gz
wget -O $1 https://zenodo.org/record/569601/files/$orig
# Need to update timestamp on file, otherwise it keeps the timestamp from
# upload day to Zenodo
touch $1

## Extract file, in the data/raw_otu_tables/ directory as well
targetdir=${target%/*}
tar -C $targetdir -xvf $target
