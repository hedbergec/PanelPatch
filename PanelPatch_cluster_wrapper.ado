program define PanelPatch_cluster_wrapper //creates k clusters and saves to generate, helps for clustersplit program below
    syntax varlist [if] [in], GENerate(namelist max = 1) k(numlist) 
    marksample touse
    capture confirm variable `generate'
    if _rc == 0 {
        drop `generate' 
    }
    
    cluster kmeans `varlist' if /// kmeans on varlist if
            `touse', ///is useable
            k(`k') measure(Gower) gen(`generate') //start(segments) //k groups, save to runningk
    
end

