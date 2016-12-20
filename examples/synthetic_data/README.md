### Synthetic data

This is a toy model of 4 atoms.

### Manifest
* `toy.pdb` - pdb file of toy system
* `toy.psf` - psf file for toy system
* `gen_psf.tcl` - tcl scripte to generate psf file with vmd psfgen
* `toy.str` - self containted str file for toy model so there's no need to load entire cgenff
* `Toy_model_MLE.ipynb` - ipython notebook of MLE of parameters for toy model.

I tried several different conditions to fit the torsion parameters of the toy model. The results are in the following
folders:
* `discrete/` - Labels (multiplicities) are on and phases discrete (0,1) with 1 representing 180.
* `discrete_decouple_n/` - Labels are off and phases are discrete (0,1) with 1 representing 180.
* `continuous/` - Labels are on and phases are continuous between 0 and 180.
* `cont_decouple_n/` - Labels are off and phases are continuous between 0 and 180.
* `neg_K_0_phase/` - Labels are on. Force constants (Ks) can take on negative values and sampling of phases is eliminated.
all phases are set to 0. A negative K value is equivalent to a phase angle of 180.
* `neg_K_decouple_n/` - Labels are off. Force constants as above.

regular MCMC has a hard time overcoming barriers when phase angle can only sample 0 or 180. The continuous phase does
a lot better. When phases are eliminated and the Ks can take on negative values the sampler also preforms well (the
mm energy is even closer to the true values than continuous phases). Labels don't seem to make a big difference. When
they are off the Ks of the multiplicities that should be off go to zero.