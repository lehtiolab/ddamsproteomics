# lehtiolab/ddamsproteomics: Configuration for other clusters

It is entirely possible to run this pipeline on other clusters, though you will need to set up your own config file so that the pipeline knows how to work with your cluster.

> If you think that there are other people using the pipeline who would benefit from your configuration (eg. other common cluster setups), please let us know. We can add a new configuration and profile which can used by specifying `-profile <name>` when running the pipeline.

If you are the only person to be running this pipeline, you can create your config file as `~/.nextflow/config` and it will be applied every time you run Nextflow. Alternatively, save the file anywhere and reference it when running the pipeline with `-c path/to/config` (see the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more).

A basic configuration comes with the pipeline, which runs by default (the `standard` config profile - see [`conf/base.config`](../conf/base.config)). This means that you only need to configure the specifics for your system and overwrite any defaults that you want to change.

## Cluster Environment
By default, pipeline uses the `local` Nextflow executor - in other words, all jobs are run in the login session. If you're using a simple server, this may be fine. If you're using a compute cluster, this is bad as all jobs will run on the head node.

To specify your cluster environment, add the following line to your config file:

```nextflow
process.executor = 'YOUR_SYSTEM_TYPE'
```

Many different cluster types are supported by Nextflow. For more information, please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html).

Note that you may need to specify cluster options, such as a project or queue. To do so, use the `clusterOptions` config option:

```nextflow
process {
  executor = 'SLURM'
  clusterOptions = '-A myproject'
}
```


## Software Requirements
To run the pipeline, several software packages are required. How you satisfy these requirements is essentially up to you and depends on your system. If possible, we _highly_ recommend using either Docker or Singularity.


### Docker
Docker is a great way to run lehtiolab/ddamsproteomics, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

Nextflow has [excellent integration](https://www.nextflow.io/docs/latest/docker.html) with Docker, and beyond installing the two tools, not much else is required - nextflow will automatically fetch the images that are needed at run time.

To add docker support to your own config file, add the following:

```nextflow
docker.enabled = true
```


### Singularity image
Unfortunately Singularity is currently not testd as a containerization method for this pipeline!
Many HPC environments are not able to run Docker due to security issues.
[Singularity](http://singularity.lbl.gov/) is a tool designed to run on such HPC systems which is very similar to Docker.

To specify singularity usage in your pipeline config file, add the following:

```nextflow
singularity.enabled = true
```

If you intend to run the pipeline offline, nextflow will not be able to automatically download the singularity images for you.
Instead, you'll have to do this yourself manually first, transfer the image file and then point to that.

First, pull the image files where you have an internet connection, then point the config file to the images:

```nextflow
singularity.enabled = true
process {
  withName youProcess {
    container = "/path/to/singularity_container.simg"
  }
}
```


### Conda
If you're not able to use Docker or Singularity, you can instead use conda to manage the software requirements. This is currently not supported, so you will have to supply the conda installs yourself.
