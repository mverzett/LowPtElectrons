import fill6371 as fill6371AOD 
import RAWBlocks
import BtoKee as BtoKeeRAW
all_samples = {
   'fill6371AOD' : fill6371AOD.input_files,
   'RAWBlocks' : RAWBlocks.input_files,
   'RAWTest' : ['/store/data/Run2017F/SingleElectron/RAW/v1/000/306/459/00000/B6410DA6-C8C5-E711-947F-02163E01A45E.root'],
   'RAWMC' : BtoKeeRAW.all_files,
   'RAWMCTest' : BtoKeeRAW.test_file
}
