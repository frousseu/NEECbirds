# NEECbirds
*Bird abundance and distribution data for emergency response.*

**Tracking data**

- Get metadata
- Build kernels 2 usage grid function

**Abundance data**

- Get metadata
- Kernels vs. rqss GAMs

##PROS and CONS of different methods

###GAMs

*PROS*

- Easily interpretable value (group size; number of birds seen when going to a site, by quantiles)

- Independent of the number of observers/checklists (but needs a certain number of observations to give a meaningful result)

*CONS*

- Can be difficult to implement

- Sensible to smoothing value

- May give extreme results when few observations are available and quantile is high (e.g. 90%)

- Need cells for predictions

###KERNELS

*PROS*

- Easier to implement and to use (e.g. risk: high, medium, low)

- Easy to output to polygons

*CONS*

- Not easy to interpret real meaning (group size, number of birds concerned?)

- Strongly depends on the amount of observers and checklist (lightly mitigated by a weigthing method based on numbers)

- Sensible to bandwith value and different parameters

##NOTES and IDEAS

- Use GAMs to predict quantrile group size in each kernel polygons

- Need to sum species by groups when doing GAMs

- Integrate ECSAS, tracking and ebird data within the same kernel output


