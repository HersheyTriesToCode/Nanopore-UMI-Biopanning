# you may have to set
# Set-ExecutionPolicy -Scope Process Unrestricted
# in your current powershell to be able to run this
docker run -it --rm --name umi-pipeline-extract-$PID -v $HOME/src:/usr/src -v $HOME/data:/usr/data umi
