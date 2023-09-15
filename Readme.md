# sdrf_convert

```
SDRF -> sdrf_convert -> CLI||config
```

## Usage
```
python -m sdrf_convert some.sdrf <tool> <additional_params_depending_on_tool>
```


## Development

### Getting started
Just enter
```
conda env create -f environment.yaml
conda activate sdrf_convert
```
and you are good to go.


### Receiving updates
Make sure to check if someone else altered `environment.yaml` or `pyproject.toml` after you pull the repository.
You might need to update
```
# conda environment
conda activate sdrf_convert
conda update -f environment.yaml
# python
pip install -e '.[dev]
```

### Need additional dependencies
* Binary dependencies (OpenMS, Comet, Yarn, NodeJS ...) will go into `environment.yaml`
* Python dependencies (pyopenms, pandas, ...) will go into `pyproject.toml`

This ensures conda does not mess up our python environment AND all python packages are in one place once we want to publish.

### Implementing a new converter
1. Add a new python file for the converter ending with `_converter.py`
2. Import the abstract converter: `from sdrf_convert import AbstractConverter`
3. Add a class definition for the new converter, inheriting the `AbstractConverter`: `class MyToolConverer(AbstractConverter):`
4. Call the super constructor in the new converters constructor
    ```py
    def __init__(self):
        super().__init__()
    ```
4. Implement the `convert()`-function which calls `init_converter()`

Example:
```py
# comet_converter.py
from sdrf_convert.abstract_converter import AbstractConverter

class CometConverter(AbstractConverter):
    # Implement the stuff used to convert the 
    # SDRF for the specific software.
    # E.g. constructor, conversions, lookups, ...
    def __init__(self):
        super().__init__()

    def convert(self, sdrf: Union[pd.DataFrame, IOBase, Path]) -> Tuple[str, str]:
        """
        Convert SDRF file to Comet CLI string and corresponding input file.

        Returns
        -------
        Any
            The converted SDRF file
        """
        self.init_converter(sdrf)
        # Do stuff
        return (cli, config)
```

### Code styling
1. Make sure to use [typed python](https://docs.python.org/3/library/typing.html). It is a bit more difficult to write but much easier to maintain and understand
2. Before submitting a pull request make sure to run pylint (`pylint sdrf_convert`) and fix your code styling