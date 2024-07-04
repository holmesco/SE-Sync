
# SE-Sync-Landmarks

This repository is an extension of the original **SE-Sync** algorithm that adds efficient computation and certification of landmark/feature locations in addition to solving [pose-graph SLAM](http://domino.informatik.uni-freiburg.de/teaching/ws11/robotics2/pdfs/ls-slam-tutorial.pdf). This is done by efficiently marginalizing over landmark measurement data, computing the certified optimal pose solution via **SE-Sync** and recovering the optimal landmark locations.

A description of the algorithm can be found in our [article](https://arxiv.org/abs/2206.12961)

The original **SE-Sync** repository can be found [here](https://github.com/david-m-rosen/SE-Sync).

## Getting Started

Currently, **SE-Sync-Landmarks** has only been implemented in MATLAB. 

Please see the **SE-Sync** repository for initial setup.

## Running SE-Sync-Landmarks

To solve for landmark locations in addition to solving pose-graph SLAM, the user must provide an additional flag to the **measurements** input structure of the original SE-Sync:

### ``measurements`` Data Structure:

| Field | Description |
|:---- | :---- |
|``edges`` | An (mx2)-dimension encoding the edges in the measurement network; ```edges(k, :) = [i,j]``` means that the $k^{th}$ measurement is of the relative transform from pose $i$ to pose/landmark $j$.  NB:  This indexing scheme requires that the states $x_i$ are numbered sequentially as $x_1, ..., x_{N}$.|
| ``R`` |  An m-dimensional cell array whose $k^{th}$ element is the rotational part of the $k^{th}$ measurement.|
| ``t`` |  An m-dimensional cell array whose $k^{th}$ element is the translational part of the $k^{th}$ measurement
| ``kappa`` | An m-dimensional cell array whose $k^{th}$ element gives the precision of the rotational part of the $k^{th}$ measurement. |
| ``tau`` |  An m-dimensional cell array whose $k^{th}$ element gives the precision of the translational part of the $k^{th}$ measurement.|
| <mark style="background:blue"> `` lmFlag `` </mark> | An mx1 dimensional array of booleans indicating whether the measurement is a landmark measurement (true) or a pose measurement (false). It is expected that the relative rotation measurement associated with landmarks is set to a 3x3 matrix of zeros.|

When defining the ``edges`` structure, the first $N_p$ states represent the poses and the last $N_m$ states represent the landmarks $(N = N_p + N_m )$. 

Once the input ``measurements`` structure is assembled as defined above, the ``SE_sync`` function can be called normally.

### Output

The output data structure is the same as in the standard **SE-Sync**. However, the output translation variables contain both pose and landmark translations in a single array. The first $N_p$ translations are associated with the poses and the last $ N_m $ translations are associated with the landmarks. 
