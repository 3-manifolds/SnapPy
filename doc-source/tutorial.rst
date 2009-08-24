========
Tutorial
========

The easiest way to learn to use SnapPy is to watch the screencasts
available on `YouTube <http://youtube.com/user/NathanDunfield>`_:

- Intro and quickstart: an 11 minute video with the basics: `Part I <http://www.youtube.com/watch?v=ezo19L-JTTI>`_ and `Part II <http://www.youtube.com/watch?v=Js4qwyIs-Oo>`_. 

- More advanced features.  (To be made)

The **key** thing to remember when using the SnapPy command shell window is
that you can explore objects using introspection and tab-completion::

     In [1]: Manifold? <hit return-key>
     ...instructions for creating a manifold...

So now we create a manifold::

   In [2]: M = Manifold("m004")

But what can we do with it?  ::

    In [3]: M.<hit tab-key>
    ...list of methods...

What does the "cover" method do? ::
     
     In [4]: M.cover? <hit return-key>
     ...description of cover method..
