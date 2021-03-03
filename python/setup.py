import setuptools

setuptools.setup(
    name="singlecell",
    version="0.0.6",
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    scripts=[
        'bin/mt_prepend_gtf.py',
        'bin/extract_name_embedded_in_gene_id.py',
        'bin/create_cellranger_ref.py',
    ],
    package_data={
        'mkref': ['create_cellranger_ref.sh'],
    },
)
