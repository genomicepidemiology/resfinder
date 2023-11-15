import yaml
import os
import sys

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '')] + sys.path
with open("src/resfinder/cge/__init__.py", 'r') as f:
    for l in f:
        if l.startswith('__version__'):
            version = l.split('=')[1].strip().strip('"')
data = {
    "package": {
        "name": "resfinder",
        "version": version
    },
    "source": {
        "url": "https://github.com/genomicepidemiology/cgecore/archive/refs/tags/{}.tar.gz".format(version),
    },
    "build": {
        "number": 0,
        "noarch": "python",
        "script": "{{ PYTHON }} -m pip install . --no-deps --ignore-installed -vvv"
    },
    "requirements": {
        "host": [
            "python",
            "kma",
            "blast"
        ],
        "run": [
            "python",
            "kma",
            "blast"
        ]
    },
    "about": {
        "home": "https://github.com/genomicepidemiology/resfinder",
        "summary": "resfinder test.",
        "license": "Apache-2.0"
    }
}

# Convert the data to YAML and print it
os.system('mkdir conda')
yaml_str = yaml.dump(data, sort_keys=False)

with open('conda/meta.yaml', 'w') as f:
    f.write(yaml_str)