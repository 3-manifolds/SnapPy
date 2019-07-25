========
Tutorial
========

The easiest way to learn to use SnapPy is to watch the screencasts
available on `YouTube <http://youtube.com/user/NathanDunfield>`_:

- Intro and quickstart: an 11 minute video with the basics: `Part I
  <http://www.youtube.com/watch?v=ezo19L-JTTI>`_ and `Part II
  <http://www.youtube.com/watch?v=Js4qwyIs-Oo>`_. 

- An hour-long demo `Practical computation with hyperbolic
  3-manifolds <http://youtu.be/j8enbAkAvdY>`_, recorded at the Thurston Memorial
  Conference.

- The SnapPy 2.0 `new feature demo <http://youtu.be/bCYe_a48viA>`_.

- `SnapPy, SageMath, and Docker are friends
  <https://icerm.brown.edu/video_archive/?play=1992>`_ with associated
  `materials (4 pages), including homework
  <https://icerm.brown.edu/materials/Abstracts/tw19-3-pods/Computation_Demonstration_2_]_Nathan_Dunfield,_University_of_Illinois,_Urbana-Champaign.pdf>`_.

- `Notes, problem sheets, and references
  <http://dunfield.info/warwick2017/>`_ from the LMS-CMI Research
  School in Warwick in 2017.

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
