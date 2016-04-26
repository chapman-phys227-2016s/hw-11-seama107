# PHYS227 HW 11

**Author:** _\<your name\>_

[![Build Status](https://travis-ci.org/chapman-phys227-2016s/hw-10-YOURNAME.svg?branch=master)](https://travis-ci.org/chapman-phys227-2016s/hw-10-YOURNAME)

**Due date:** 2016/05/03

## Specification

Complete the following exercises, placing your solutions into separate files. In each file, write the solution as a callable function, so that you can write suitable test functions that demonstrate correct output using the nose framework. GitHub will automatically run your tests on every commit, indicating any failures via the Travis framework with build status above.

1. Complete CW11 and use its module ```sombrero.py```.
   
1. For $F = 0.4$, create a high-resolution scatterplot of $(x(t),y(t))$ for the specific points $t = (\phi + n)2\pi$ for $\phi = 0$ and $n = 0,1,\ldots,1000$. Label your horizontal axis $x$ and your vertical axis $x'$. Plot your axes on a fixed scale $x\in[-1.5,1.5]$ and $y\in[-1.1,1.1]$.  Export the plot directly to a file ```frame00.png```.
   
1. Repeat the previous problem, but using $\phi = 0.01$, saving the file in ```frame01.png```. Repeat until $\phi = 0.99$ and ```frame99.png```. 
   
1. Create an animated gif out of the previous 100 frames, using the ```convert``` commandline tool, as before.
   
Finally, to cleanly present your work, create a Jupyter notebook ```hw11.ipynb``` that loads your animated gif and explains what sort of phase-plane motion you observe from the dynamics as time evolves. The notebook should be formatted professionally using Markdown headings, including your name and a brief description of each homework section.

## Assessment

_\<Analyze what you found useful about this assignment in your own words. Include any lingering questions or comments that you may have.\>_

## Honor Pledge

I pledge that all the work in this repository is my own with only the following exceptions:

* Content of starter files supplied by the instructor;
* Code borrowed from another source, documented with correct attribution in the code and summarized here.

Signed,

_\<your name\>_
