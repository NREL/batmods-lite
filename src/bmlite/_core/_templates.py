import os

import bmlite as bm


def templates(model: str, file: str | int = None) -> None:
    """
    Print simulation templates.

    If no input, a list of available templates will print. Otherwise, given
    a file name or index, a template will print.

    Parameters
    ----------
    model : str
        Model package name, e.g., 'SPM'. Input is case insensitive. All model
        names are converted to uppercase within this function.
    file : str | int, optional
        File name or index. The default is None.

    Returns
    -------
    None.

    Raises
    ------
    AttributeError
        'model' is not a valid model package.
    FileNotFoundError
        'model' has no 'templates' directory.

    """

    try:
        path = getattr(bm, model.upper()).__path__[0]
    except AttributeError:
        raise AttributeError(f"{model=} is not a valid model package.")

    if not os.path.exists(path + '/templates/'):  # pragma: no cover
        raise FileNotFoundError(f"{model=} has no 'templates' directory.")

    templates = os.listdir(path + '/templates/')

    if file is None:

        print('='*30, model + ' templates:', '='*30, sep='\n')
        for i, f in enumerate(templates):
            print('  - [' + str(i) + '] ' + f.removesuffix('.yaml'))

    elif isinstance(file, str):
        file = file if '.yaml' in file else file + '.yaml'

    elif isinstance(file, int):
        file = templates[file]

    if file is not None:

        print('='*30, file, '='*30, sep='\n')
        with open(path + '/templates/' + file, 'r') as f:
            print(f.read())
