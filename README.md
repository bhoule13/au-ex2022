## Effect of the COVID-19 pandemic on life expectancy in Australia

### Introduction
------------

This GitHub repository contains the replication code and data to accompany *Effect of the COVID-19 pandemic on life expectancy in Australia*.

### Prerequisites
------------

To run the code locally you will need a working installation of [**R**](https://www.r-project.org/), including the necessary dependencies.

Databases from the Australian Bureau of Statistics (ABS) are used in the analysis and included in six input data sets:

* _PopO.csv_ Estimated Resident Population (ERP) by age, sex and Australian state/territory of residence for years 2010 to 2021. The data is presented exactly as downloaded from the [ABS tabulator](https://explore.data.abs.gov.au/).
* _PopQ.csv_ Estimated Resident Population (ERP) by age, sex and Australian state/territory of residence for the year of 2022. The data is presented exactly as downloaded from the [ABS tabulator](https://explore.data.abs.gov.au/).
* _DeathsO.csv_ includes the number of deaths by age, sex and Australian state/territory of residence by year of occurrence, for years 2017 to 2021. The data is presented exactly as downloaded from the [ABS tabulator](https://explore.data.abs.gov.au/).
* _DeathsStates2.csv_ includes the number of deaths by age-group, sex and state/territory of residence occurring in the year 2022. The data were purchased from the ABS.
* _CausesNational3.csv_ includes the proportion of causes of death by 5-year age group, sex for Australia. The selected causes of death include: COVID-19, Influenza, Pneumonia, Other respiratory diseases, Neoplasms, Ischaemic heart disease, Cerebrovascular diseases, Dementia including Alzheimers, Diabetes, External causes, Other causes, and External causes and Other causes. Data from the World Health Organization (WHO) Mortality Database was used for the data by cause and age in the years 2017-2020. The data for 2021 and 2022 were purchased from the ABS
* _CovidState2.csv_ includes the proportion of deaths by age, sex, state/territory of residence, and COVID-19 and other causes of death, occurring in the year 2022. The data were purchased from the ABS.

Death rates and population under exposure by age, and life expectancy at birth for females and males from 27 countries from the Human Mortality Database (HMD: [www.mortality.org](www.mortality.org)) are used. We recommend users downloading the most recent information directly from this data source.

### Structure
----------------

* _r_ contains R code. The code assumes the working directory is set to the _r_ directory.
* _data_ contains the data for Australia.
* _output_ contains outputs (i.e., figures) produced from the R code.

### Acknowledgements

The R code and resulting figures partly re-use [publicly available code and data](https://zenodo.org/record/6861804) used by *Schöley et al. (2022)*, with some modifications.

Schöley, Jonas, Kashnitsky, Ilya, Kniffka, Maxi S., Zhang, Luyin, & Kashyap, Ridhi. (2022). Analysis code for "Bounce backs amid continued losses: Life expectancy changes since COVID-19" (v0.98). Zenodo. https://doi.org/10.5281/zenodo.6861804