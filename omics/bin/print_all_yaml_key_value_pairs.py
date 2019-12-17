#!/usr/bin/env python3

import yaml
import sys


if len(sys.argv) != 2:
    print("Usage: " + sys.argv[0].split("/")[-1] + " YAML_FILE")
    sys.exit(1)


def get_all_key_value_strings(concatenated_key, element):
    key_base = concatenated_key
    if isinstance(element, dict):
        for key, value in element.items():
            concatenated_key = key_base + "_" + key
            get_all_key_value_strings(concatenated_key, value)
    elif isinstance(element, list):
        for item in element:
            get_all_key_value_strings(concatenated_key, item)
    else:
        print(concatenated_key + "=" + str(element))
    return


data = yaml.load(open(sys.argv[1], "r"), Loader=yaml.FullLoader)
for key, element in data.items():
    concatenated_key = key
    get_all_key_value_strings(concatenated_key, element)
