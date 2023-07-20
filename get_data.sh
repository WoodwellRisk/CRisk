# Make STORM directory in input
mkdir input/STORM

# Copy STORM track files from bucket to local environment
gsutil -m cp gs://cmip5_data/tropical_cyclones/STORM/tracks/*$1*.txt ./input/STORM