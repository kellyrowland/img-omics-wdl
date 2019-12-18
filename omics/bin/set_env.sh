#!/usr/bin/env bash

if [[ $# -ne 1 ]]
then
	echo "Usage: $(basename $0) <run_config.yaml>\n" >&2
	exit 1
fi

config_yaml=$1
if [[ ! -f $config_yaml ]]
then
	echo "$config_yaml is not a file! Aborting!" >&2
	exit 1
fi

echo "$(date +%F_%T) - Creating env variabes from yaml conf now..."

set -a   # Sets the export attribute for each variable between -a and +a.
eval $(./print_all_yaml_key_value_pairs.py $config_yaml)
# Create specific /usr/bin/time output format.
TIME="\n-----/usr/bin/time output-----\n"
TIME=$TIME"Cmd: %C\nCPU: %P\nMemory: %M KB (res max)\n"
TIME=$TIME"Time: %e (%E real, %U user, %S sys)\nExit code: %x\n"
set +a

echo "$(date +%F_%T) - Done exporting the environment variables."
