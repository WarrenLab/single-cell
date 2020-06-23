import setuptools

setuptools.setup(
    name="singlecell",
    packages=setuptools.find_packages(),
    install_requires=['python>=3.6'],
    scripts=[
        'mt_prepend_gtf.py',
        'create_cellranger_ref.py',
    ]
)
