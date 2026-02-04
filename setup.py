from setuptools import setup, find_packages

setup(
    name="mutme",
    url="https://github.com/vidrl/mutme",
    author="Eike J. Steinig",
    author_email="eike.steinig@mh.org.au",
    packages=find_packages(),
    # include_package_data=True, 
    # package_data={
    #     'utils.assets': ['ercc.tsv']
    # },
    install_requires=[
        "typer",
    ],
    entry_points="""
        [console_scripts]
        mutme=mutme.terminal:app
    """,
    version="0.1.0",
    license="MIT",
    description="Generic mutation and mutation-of-interest calling from consensus genomes",
)