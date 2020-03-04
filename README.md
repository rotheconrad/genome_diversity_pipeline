# genome diversity pipeline

To add a dataset:

```bash

/path/to/pkg/01_add.bash target_dir dataset_name /path/to/reads.[12].fastq.gz

```

This will automatically copy and launch the pipeline step by step.


If you want to continue a halted run:

```bash

/path/to/pkg/xx_runme.bash target_dir dataset_name

```

This will detect the last complete step and lanch the next one.

