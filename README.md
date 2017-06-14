# Event Lifetime

The lifetime of an event is the time that it takes for the moving brightness gradient causing the event to travel a distance of 1 pixel.
The provided algorithm augments each event with its lifetime, which is computed from the event's velocity on the image plane.
The generated stream of augmented events gives a continuous representation of events in time, hence enabling the design of new algorithms that outperform those based on the accumulation of events over fixed, artificially-chosen time intervals.
A direct application of this augmented stream is the construction of sharp gradient (edge-like) images at any time instant.

![building](https://user-images.githubusercontent.com/670994/27136305-48c4b81c-511b-11e7-8256-ac7288b24783.png)

For more details, please read our [ICRA'15 paper](http://rpg.ifi.uzh.ch/docs/ICRA15_Mueggler.pdf).

## Disclaimer and License

This code has been tested with MATLAB R2016b on Ubuntu 16.04.
This is research code, expect that it changes often and any fitness for a particular purpose is disclaimed.
The source code is released under a GNU General Public License (GPL).


## Instructions

Please run the file `matlab/main.m`.

It computes the lifetime of each events and creates a video with sharp renderings and, as comparison, another video using a fixed integration time.
Sample data is included in this repository.


## Publication

If you use this code in an academic context, please cite the following [ICRA'15 publication](http://rpg.ifi.uzh.ch/docs/ICRA15_Mueggler.pdf):

E. Mueggler, C. Forster, N. Baumli, G. Gallego, D. Scaramuzza:
**Lifetime Estimation of Events from Dynamic Vision Sensors.**
IEEE International Conference on Robotics and Automation (ICRA), Seattle, 2015.

    @inproceedings{Mueggler15ICRA,
      author = {Mueggler, Elias and Forster, Christian and Baumli, Nathan and Gallego, Guillermo and Scaramuzza, Davide},
      title = {Lifetime Estimation of Events from Dynamic Vision Sensors},
      booktitle = {IEEE International Conference on Robotics and Automation (ICRA)},
      year = {2015}
    }
