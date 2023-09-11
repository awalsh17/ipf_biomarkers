These files are called from the {targets} pipeline 

+ functions.R
+ main_report.Rmd

```
targets::tar_make()
```

These files need to be run first. Yes, they should have just been added to the top of the pipeline.

+ download_geo_data.R
+ sc_analysis.R
