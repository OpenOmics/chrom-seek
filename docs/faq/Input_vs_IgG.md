# Input control versus IgG

IgG controls and input controls in ChIP are two very different types of controls and 
cannot be used interchangeably in most projects. Input controls are where you sonicate 
chromatin and directly turn that into a library so it accounts for openness of the chromatin. 
This is the kind of control that MACS can correct for directly in its model. IgG controls are 
pull-downs with IgG agarose beads without any antibodies. These libraries tend to be sequenced at 
much lower depth than all other samples, are typically super sparse, typically break all QC metrics 
due to over PCR-amplification, and in many cases are not overall useful.
