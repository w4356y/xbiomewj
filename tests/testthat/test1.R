context("Test Plot Function")

test_that("Test my foo function: ",
      {
            expect_error(plot_pcoa(otu_table(phylo), sample_data(phylo), "Response", "bray"), NA)
            expect_error(plot_pca(otu_table(phylo),sample_data(phylo),"Response"), NA)
            expect_s3_class(permanova_test(otu_table(phylo), sample_data(phylo), "Response","bray"), "adonis")
            expect_error(advanced_permanova_test(otu_table(phylo), sample_data(phylo), "Response","Diagnosis","bray"), NA)
            expect_error(plot_pcoa_with_shape(otu_table(phylo),sample_data(phylo),"Response","Diagnosis","bray"), NA)
            expect_error(plot_pcoa_plus_seq(plot_pcoa(otu_table(phylo), sample_data(phylo), "Response", "bray"),otu_table(phylo) %>% as.data.frame(), Response), NA)
            expect_error(plot_nmds(otu_table(phylo), sample_data(phylo), "Response"), NA)
      })
