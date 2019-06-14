These are exemplary steering files to calculate NLO predictions for the total and single-differential 
ttbar cross sections using MNR calculatinos intefraced in xFitter. They work with PionCeres branch 
of xFitter as of 22.05.2019. 

Assume that you have just installed xFitter PionCeres. Create a new directory e.g. 'run' and unpack 
input files there:
mkdir run
cd run
tar xvzpf mnr-steering.tar.gz

Look at the unpacked files. The main steering files which xFitter needs for any run are steering.txt, 
parameters.yaml and constants.yaml. They might be very complex for the first look. However, try to compare 
with those in the main xfitter directory, perhaps it is 'xfitter-master', run e.g. 
'diff xfitter-master/steering.txt run/steering.txt'
and see the changes. All the changes are related to the MNR calculatinos, and you need to understand 
them in order to be able to change anyhting in these calculations in the future, while for the time being you may 
ignore all other settings in the input files which are identical in xfitter-master and mnr-steering.tar.gz. 

This is a brief summary of the changes:

steering.txt
  - only three input data files are enabled (for the total, pT and y cross sections); the data files will be explained later
  - disabled writing LHAPDF6 grids (no need to understand, just a technicality to save disk space)

parameters.yaml
  - PDFs are taken from 'ABMP16_5_nlo' PDF set (b.t.w. you should have it installed)
  - settings for the 'cbdiff' reaction (this is the one which runs MNR) are provided:
    -- Order: NLO (should be clear, could be LO or NLO)
    -- MS_MASS: 0 you use 0 for pole mass or 1 for MSbar mass
    -- mr_A: 1.0, mr_B: 1.0 this is how the renormalisation scale is calculated:
       mu_ren = sqrt(mr_A * p_T^2 + mr_B * m_Q^2) where p_T is transverse momentum of the heavy quark
       and m_Q is its mass, i.e. we have mu_ren = sqrt(p_T^2 + m_Q^2)
    -- mf_A: 1.0, mf_B: 1.0 same but for the factorisation scale
    
constants.yaml
  - no changes :)

Now take a look at the data files mnr-tt-total.dat, mnr-tt-pt.dat, mnr-tt-y.dat. 
The user specifies in such files what should be calculated. Let's look at mnr-tt-total.dat:
  - Name = "Dummy", IndexDataset = 999 are irrelevant now
  - TermName, TermSource and TermInfo are arrays of theory "terms" (each term can be a separate process).
    Now we have one term 'A' (the name is arbitrary) which is the 'cbdiff' reaction (MNR). TermInfo supplies 
    steering and is the most important part. The TermInfo settings are in format key=value separated by :
    Here they are described:
     - steer_q=0:steer_a=1 enable predictions for antiquarks only (we needed this to compare to MCFM where by 
       default antitop distributions wee calculated, though the difference betweed q and aq distributions is tiny)
     - mq=173.0 heavy quark mass in GeV
     - FinalState=parton:FragFrac=1.0:FragPar=1.0 should be always like this so far (could be c or b hadrons with 
       the extra fragmentation function and fraction applied, but so far the final state is partons)
     - energy=8000.0 c.m.s energy in GeV
     - pTn=1:pTmin=0.0:pTmax=1500.0:yn=1:ymin=0.0:ymax=4.0 the reaction calculates double-differential cross sections
       as function of pT and y, and here you specify the bins. For the total cross section, I specified one bin 
       for both pT and y with sufficiently large upward boundaries. Look at mnr-tt-y.dat, mnr-tt-pt.dat to see how 
       one specifies bins for single-differential cross sections.
     - steer_ptmin=0.001:steer_ptmax=1600.0:steer_npt=50:steer_nptsm=100 now the tricky settings come. This is an internal
       pT binning (the names have 'steer_' prefix) used when calculating integrals. npt is the number of bins, typically should be 
       ~ 100, min is the lowest boundary which should be very small (but not 0 to avoid numerical problems),
       max is the highest boundary (should be a bit larger than pTmax, because there is cubic spline interpolation applied), 
       nptsm is the number of bins after smooth interpolation (should be ~2 more than npt)
     - steer_ymin=0.0:steer_ymax=4.0:steer_ny=50 same for y, but there is no smooth interpolation applied when calculating integrals
     - steer_nsfnb=500 number of bins of PDFs as function of partonic x (also called "structure functions"), 
       results should be barely sensitive to this number
     - steer_nx3=25 number of bins as function of the rapidity boost of the parton-parton system, should be optimisaed for different 
       heavy quakr masses, kinematics etc.
     - steer_nx4=125 number of bins in variable which controls NLO radiation (relevant for NLO corrections value, 125 works 
       reasonably well for all calculations I've done)
     - steer_nbz=50 relevant only if you specified FinalState!=parton and apply a fragmentation function
     - steer_xmin=1e-3:steer_xmax=1.0 the range of partonic momentum fraction x (depends on the heavy quark mass and kinematics
       of the final state)
     - steer_mf2min=2500.0:steer_mf2max=11240000.0 the range of factorisation scale squared [in GeV^2], depends on the settings 
       mf_A, mf_B, heavy-quark mass and pT interval (a good idea is to choose the minimum and maximum values four times lower 
       and upper then estimated for the nominal scale choice, because later you might want to vary the scale up and down by factor 
       2; in general the results are barely sensitive to the interval of mf2min, mf2max even if it is too large)
     - steer_nf=5 number of light flavours; MNR calculations treat quarks either as massless or massive, e.g. for top production 
       one treats u, d, s, c and b as massless; this setting should be consistent with the chosen PDF set, therefore we have 
       'ABMP16_5_nlo' with 5 flavours; be aware that you need 3 for charm and 3 or 4 for beauty (and corresponding PDF sets).
     - furthermore, any parameter from parameters.yaml cbdiff reaction, like MS_MASS etc., can be overriden for individual 
       data sets by specifying it in TermInfo. E.g. if for whatevet reason one wants to change the order to LO for the y cross section
       only, one adds Order=LO in TermInfo in mnr-tt-y.dat (the other files will keep using Order NLO specified in parameters.yaml)
  - TheorExpr = 'A*2' this is the final expression for theory prediction. In this case the MNR cross section will be in microbarn 
    times 2 because we set ymin=0, ymax=4, while we need the cross section from -4<y<4 (essentially -inf<y<inf), and in pp collisions
    the cross section is symmetric for y<0 and y>0, therefore no need to calculate y<0
  - NData =   1 the number of "data" points. One needs to supply "data" (which are rather pseudodata) because xFitter is a fitting framework 
    which always calculate theory predictions and compare to data in order to extract theory parameters. You will use only the calculated 
    theory predictions so the fake data are irrelevant.
  - NColumn =   4
    ColumnType = 2*"Bin","Sigma", 1*"Error"
    ColumnName = "Bin_min","Bin_max","Sigma", "stat"
    Percent = 1*true
    these settings you do not change, they are needed to read the fake data correctly (two bin boundaries, cross section and statistical
    uncertainty in percent)
  - &PlotDesc..&End these is steering for the plotter, One uses any bin as x-axis, PlotN=NData, the title of the axes are specified
  - and then the fake data values come. The only thing is the number of lines should match NData (look at mnr-tt-pt.dat, mnr-tt-y.dat 
    where more data points have to be supplied)

To run the code, type
xfitter

It should take less then one minute to calculate the cross sections. The output is in output/fittedresults.txt (the name is misleading, 
this is not a fit of course). Compare this file to the reference one which I produced [fittedresults.txt.reference]. 
The calculated theoretical predictions for each bin are in the 7th column. If double-differential pT-y cross sections are requested, 
the order is: all pT bins in the 1st y bin, then all pT bins in the 2nd y bin etc.

If you want to look at the calculations themselves, the relevant files (with a lot of code) are
  - reactions/cbdiff/src/Reactioncbdiff.cc
  - reactions/BaseHVQMNR/src/ReactionBaseHVQMNR.cc
  - reactions/BaseHVQMNR/src/hvqcrsx.f
  - reactions/BaseHVQMNR/src/MNR.cc
  - reactions/BaseHVQMNR/src/MNRFrag.cc
  - reactions/BaseHVQMNR/src/MNRGrid.cc

Now you need to make sure that all steer_... settings in TermInfo were fine, and the precision of calculated cross sections is reasonable.
For this the way I am using is varying them up and down one setting at a time, re-running the calculation and looking at the differences 
in each bin. If I find the differences are significant (I would say <0.1% is not significant), then I change the setting and 
repeat the procedure (yes, it is tedious a bit). The guidelines are:
 - with too many bins the program will be too slow and will use too much memory (look at the messages like 'MNR::Precalc(): required 19 MB' 
   for each data set in the console when running xfitter, you certainly do not want >1GB memory usage)
 - with too few bins results will be imprecise
 - with too wide ranges (e.g. for internal pT, y integrals, or for mf2 and x) results will be imprecise, because e.g. the calcultion will 
   be done at very high pT and y where the cross section is negligibly small
 - with too narrow ranges, e.g. if you set steer_ymax=2.5, for the total cross section the contribution at y>2.5 is not calculated, 
   while it could be a few % of the total cross section
so you need to find some balance between the number of bins and ranges on one hand, and required accuracy and running time/memory usage on 
other hand. The numbers in these steering should be already optimised, but you can train a bit on them to get an idea how to do it for 
another calculation in the future. As an example, try to vary steer_nx4 by changing it to e.g. 50 or 250 and see how the running time 
and memory usage are changing (some settings do not affect memory usage), and then see how the results are changing. Before running with 
the changed setting, specify another output directory by adjusting OutDirName = 'output' in steering.txt (to not override the old results). 
Furthermore, there is a handy way to monitor the changes in the cross sections by plotting the results (assumed you changed nx4 to 50 and 
specified OutDirName = 'output-nx4-50'):
xfitter-draw --no-pdfs --no-shifts --no-tables --only-theory --outdir plots-nx4-50 output output-nx4-50
Look at plots-nx4-50/plots.pdf to see the cross sections and their ratios (the refernce plot which you should obtain is 
plots-reference.pdf). You could see the changes are about 0.4%. This I would call significant, so I do not want to change steer_nx4 to 50.
But maybe 100 would be fine? Check yourself! Then try increasing it. E.g. you can scan (just arbitrary numbers):
250 vs 200 difference ~0.01%
250 vs 150 difference ~0.02%
250 vs 125 difference ~0.05%
250 vs 100 differernce ~0.1%
250 vs 50 difference ~0.4%
Clearly nothing changes if you increase it above 250. For such numbers, 125 looks fine to me. But in practice I keep even higher values 
if running time < 30s and memory usage < 100MB.
This game might be tricky if there are differntial cross section with many bins. You will need to decide whether you really want 0.1% 
accuracy in all bins (e.g. for bins with very low cross section such as 3<y<4 it does not make sense). This was a long part, but in the 
future when you will produce calculations for c and b quarks, you will need to changes these steering numbers, and you need to know how 
to do it, and this is the most important part of the game: to make sure you use reasonable settings in TermInfo.
