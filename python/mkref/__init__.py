"""
Functions and data for running cellranger mkref
"""
import pkg_resources


def make_sbatch_command(sbatch_params: dict, sbatch_exports: dict):
    """
    Constructs an sbatch command to mkref

    Given sbatch parameters and sbatch exports, finds the resource
    'create_cellranger_ref.sh' and constructs an sbatch command,
    which is then returned as a list of strings.

    Args:
        sbatch_params: dict whose keys are the names of sbatch
            flags and values are the values to be passed to those
            command-line arguments. E.g., the pair ('mem', '10G')
            becomes '--mem=10G'.
        sbatch_exports: dict whose keys are the names of variables
            to be exported by sbatch to the bash script, and values
            are the values to set those exported variables to

    Returns: a list of strings containing a command that can be used
        to run sbatch, compatible with methods such as subprocess.run
    """
    script_path = pkg_resources.resource_filename(
        __name__, 'create_cellranger_ref.sh'
    )
    params_list = [f'--{k}={v}' for k, v in sbatch_params.items()]
    exports_string = ','.join([f'{k}={v}' for k, v in sbatch_exports.items()])
    params_list.append(f"--export='{exports_string}'")
    return ['sbatch'] + params_list + [script_path]
