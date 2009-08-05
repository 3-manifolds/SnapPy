========
Tutorial
========

The easiest way to learn to use SnapPy is to watch the screencasts
available here:

- Intro and quickstart. (To be made)
- More advanced features.  (To be made)

The **key** thing to remember when using the SnapPy command shell window is
that you can explore objects using introspection and tab-completion::

     In [1]: Manifold? <hit return-key>
     ...instructions for creating a manifold...

So now we create a manifold::

   In [2]: M = Manifold("m004")

But what can we do with it?  ::

    In [3]: M.<hit tab-key>
    42 possibilities 
    <hit tab-key again>
    ...list of methods...

What does the "cover" method do? ::
     
     In [4]: M.cover? <hit return-key>
     ...description of cover method..
