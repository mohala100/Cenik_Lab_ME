# https://ribosomeprofiling.github.io/ribor/ribor.html
# https://watermark.silverchair.com/btaa028.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAsAwggK8BgkqhkiG9w0BBwagggKtMIICqQIBADCCAqIGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMdbGgSSuD7arvaO6GAgEQgIICcywgWDJBd6HiD1oBJoqlZ_KbxhwRt9Z73kCXxvitrOcHVMbprB0-K9vX9ExmHbb0-2CktuR03sKjExE33u3RN81MTMTDoMw2itfz5tLzDe6fRpnPd2RTY51gszmEAQYwqb7mDxGE52J49WdUBoyQ5qkQiTG_RGdiJPtbq5Cs7nen2C8xrjusA1F8V4CK8HsSThbO8XD0eNb6PMU1rCyobMVvCzPDCxg9gg6DpNWg4XY4h8Ggt7noxZopYwwXZThhEqqSr4RQgjZ2tuP5QoPRlkQnH7Es8dnouWjwQ7t304YKjdgow1Ikwl6wYoV4jn4wQocBpHfJkeAqG8bmOrvzjBZ0f_HwWqvZNo-SakfyVxwu0Rm12fu9xKVmr22bzFnwIM0-evRChFCxbxEvUPjnemnHaFCyqLFqG6wRE_FKKW8d1ySz2oUOLC-bIEO5ouLhXWx2LHv3kUtIvaistC_t0fk0uqI8FdOLfBpA_9M-lxDFJwRJmf73eM_MckcyX6xCop40A40e2IWSV-HvYlfjmywypOdwpyjhjCNRIi0KxsqUFJgdCcHLkj1Q0-1Q1fd5NZ6UvXAX0CPdnYUOMG4B-xKFO1e780CimjrPLnfth3VLDWfGnvwCjNC0LZjbdxYxISV3NWtye6Zs0xYmkQHJJ0rxhGl3QCZXZlqnil2M_94vBFLTHC_ZAJawM9XfT-kvRPxlMqK605qbpO6CvhGH9SoX6UUZLzwH3Z-XyS9EiPY0E9qOfc7cNcKYJeYkS0ENnwRvRlfOs4Q04J0MjU0is40f2B36hcTrh4dU_BhD_MseToSOTHGd3ImFSZS_LsfafgqxUg
# https://ribosomeprofiling.github.io/
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ribor")

install.packages("devtools")
library("devtools")
install_github("ribosomeprofiling/ribor")

library(ribor)


#file path to the example ribo file 
file.path <- system.file("extdata", "HEK293_ingolia.ribo", package = "ribor")

#generates the 'ribo' class object 
original.ribo <- Ribo(file.path, rename = rename_default )
original.ribo

plot_length_distribution(x           = original.ribo,
                         region      = "CDS",
                         range.lower = 28,
                         range.upper = 32,
                         fraction    = TRUE)

plot_length_distribution(x           = original.ribo,
                         region      = "CDS",
                         range.lower = 28,
                         range.upper = 32,
                         fraction    = FALSE)

rc <- get_length_distribution(ribo.object      = original.ribo,
                              region      = "CDS",
                              range.lower = 28,
                              range.upper = 32)
rc

plot_length_distribution(rc, fraction = TRUE)

get_info(original.ribo)$attributes$metagene_radius

plot_metagene(original.ribo,
              site        = "start",
              experiment  = c("GSM1606107"),
              range.lower = 28,
              range.upper = 32)

plot_metagene(original.ribo,
              site        = "stop",
              normalize   = TRUE,
              title       = "Stop Site Coverage",
              range.lower = 28,
              range.upper = 32)

meta.start <- get_metagene(ribo.object = original.ribo, 
                           site        = "start",
                           range.lower = 28,
                           range.upper = 32,
                           length      = TRUE,
                           transcript  = TRUE)

print(meta.start[ , 1:10])

tidy.meta.start <- get_tidy_metagene(ribo.object = original.ribo,
                                     site        = "start",
                                     range.lower = 28,
                                     range.upper = 32,
                                     length      = TRUE)
head(tidy.meta.start, 2)
