venv_cmd=mamba

# Create environment
$venv_cmd create -n $1 python=3.9
$venv_cmd activate $1

# Install CLIMADA
$venv_cmd env update --file env_surge.yml
pip install -e .