# you may have to set
# Set-ExecutionPolicy -Scope Process Unrestricted
# in your current powershell to be able to run this
docker run --platform linux/x86_64 -it --rm --name pipeline-umi-amplicon-$PID -v ~/src:/usr/src -v ~/data:/usr/data pipeline-umi-amplicon
