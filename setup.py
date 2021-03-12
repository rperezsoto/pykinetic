from setuptools import setup, find_packages

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name = 'pykinetic2',
    packages = ['pykinetic2',],
    description = """A python library and command line apps to write 
                     microkinetic models for """,
    keywords = ['compchem', 'microkinetics'],
    author = 'Raúl Pérez-Soto',
    classifiers = ["Programming Language :: Python :: 3",],
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url = 'https://github.com/fmaserasgroup-repo/pyssian-utils', # Update
    python_requires='>=3.6',
    install_requires=['setuptools','pathlib','numpy'],
    extras_requires=['numpy','scipy'],
    include_package_data=True,
    package_data = {'templates': ['pykinetic2/templates/*'],
                    'tests'    : ['pykinetic2/tests/*.py']},
    scripts = ['pykinetic2/scripts/pykinetic-model.py',
               'pykinetic2/scripts/pykinetic-scan.py'],
)
