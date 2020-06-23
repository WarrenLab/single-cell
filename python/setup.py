import setuptools

setuptools.setup(
    name="singlecell",
    packages=setuptools.find_packages(),
    scripts=[
        'mt_prepend_gtf.py',
        'create_cellranger_ref.py',
    ]
)
