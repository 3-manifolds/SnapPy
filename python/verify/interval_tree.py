from __future__ import print_function

__all__ = ['IntervalTree']

LEFT  = 0
RIGHT = 1

class IntervalTree(object):
    """
    A data structure that can store (interval, value) pairs and quickly
    retrieve all values for which the interval overlaps with a given
    interval.

    Create an interval tree and add pairs to it::

        sage: from sage.all import RIF
        sage: t = IntervalTree()
        sage: t.insert(RIF(1.01,1.02),'1')
        sage: t.insert(RIF(3.01,3.02),'3')
        sage: t.insert(RIF(2.01,2.02),'2')
        sage: t.insert(RIF(0.99, 3.0),'big')

    Retrieve all values for intervals overlapping [1.5, 2.5]::

        sage: t.find(RIF(1.5, 2.5))
        ['big', '2']

    Calls to insert :py:meth:`insert` and :py:meth:`find` can be mixed::

        sage: t.insert(RIF(4.01,4.02),'4')
        sage: t.find(RIF(1,10))
        ['big', '1', '2', '3', '4']
        sage: t.find(RIF(16, 17))
        []

    This is implemented as interval tree, i.e., as a red-black tree keyed by
    the left endpoint of an interval with each node *N* storing the max of all
    right endpoints of all intervals of the subtree under *N*.

    A demo of red-black trees is at https://www.youtube.com/watch?v=gme8e_6Fnug
    and explanation at http://unhyperbolic.org/rbtree.pdf .
    Also see wikipedia or Cormen et al, *Introduction to Algorithms*.
    """

    class _Node(object):
        def __init__(self, interval, value):
            # The interval and value of the inserted pair this node corresponds
            # to. The left endpoint of the interval is used as key for the
            # binary search tree.
            self.interval = interval
            self.value = value
            # The maximum of all right endpoints of all intervals of the
            # subtree under this node, including the node itself.
            self.max_value = interval.upper()
            # Left and right child
            self.children = [None, None]
            # Is this node red or black?
            self.isRed = True

        def update_max_value(self):
            # Assuming that the children store the correct max_value, update
            # this node's max_value
            self.max_value = self.interval.upper()
            for child in self.children:
                if child:
                    self.max_value = max(self.max_value, child.max_value)

    def __init__(self):
        self._root = None

    def find(self, interval):
        """
        Finds all values that have been inserted with intervals overlapping the
        given interval which is an element in SageMath's ``RealIntervalField``.

        The runtime of this call is O(log (n) + m) where n is the number of
        intervals in the tree and m the number of the returned intervals.

        The order of the returned values is ascending in the left endpoint of
        the associated intervals. If several intervals have the same left
        endpoint, the order of the returned values depends on the insertion
        order but is still deterministic.
        """

        result = []
        # Recurse the tree
        IntervalTree._fill_recursive(self._root, interval, result)
        return result

    @staticmethod
    def _fill_recursive(node, interval, result):
        if not node:
            return

        # The max_value is the maximum over the entire subtree. If it is
        # strictly lower than the queried interval, we can stop recursing.
        if node.max_value < interval.lower():
            return

        # Recurse down to the left
        IntervalTree._fill_recursive(node.children[LEFT], interval, result)

        # The left endpoint of the interval is used as key. If it is strictly
        # greater than the queried interval, we do not need to consider this
        # node or the subtree to the right.
        if node.interval.lower() > interval.upper():
            return

        # Add to result if overlapping
        if node.interval.overlaps(interval):
            result.append(node.value)

        # Recurse down to the right
        IntervalTree._fill_recursive(node.children[RIGHT], interval, result)
        
    def insert(self, interval, value):
        """
        Inserts (interval, value) where interval is an element in
        SageMath's ``RealIntervalField`` and value can be anything.
        """

        # Create node
        node = IntervalTree._Node(interval, value)

        if self._root:
            # Recurse down the tree to insert node, fix the red-black
            # property and update max field when unwinding the stack.
            #
            # The root of the tree might change, so update self._root.
            violation, self._root = IntervalTree._insert_fix_and_update_max(
                self._root, node)
        else:
            # If empty tree, node becomes root
            self._root = node
        
        # Root is always black
        self._root.isRed = False

    @staticmethod
    def _if_red(node):
        # Does node exist and is red
        if node and node.isRed:
            return node
        return None

    @staticmethod
    def _insert_fix_and_update_max(node, leaf):
        # We actually recurse down the tree to figure out where to insert a
        # node and then use the unwinding of the stack to fix the red-black
        # tree and update the max_value field.
        #
        # This functions does this in three steps:
        # Step 1: figure out whether to follow the left or right child and
        #         call yourself to recurse down the tree to insert the new
        #         node.
        # Step 2: fix any potential violation where a black node has a red
        #         child and a red grand child.
        # Step 3: update the max_value field of the node.
        #
        # The result of this function is a pair (violation, node) used by step
        # 2 of the calling function is:
        # * violation is
        #       - 'VIOLATION' if there is definitively a violation where a red
        #         node has a red child
        #       - 'POTENTIAL' indicating that the node has became red and there
        #         is a violation if and only if the parent is red
        #       - 'OK' meaning there is no violation
        # * node is the new root of the subtree - the root of the subtree might
        #   have changed due to rotations so we need to update the respective
        #   child field of the parent to point to the new root.

        # Step 1 and 2
        violation, node = IntervalTree._insert_and_fix(node, leaf)

        # Step 3
        node.update_max_value()
        
        return violation, node

    @staticmethod
    def _insert_and_fix(node, leaf):
        # Step 1: find whether to follow left or right child
        
        branch = LEFT if leaf.interval.lower() <= node.interval.lower() else RIGHT

        # If we recursed down and have no more child, ...
        if not node.children[branch]:
            # ... add this as leaf node. The new leaf node is always red, ...
            node.children[branch] = leaf
            # ... so there is something to fix if and only if this node is also
            # red.
            return ('VIOLATION' if node.isRed else 'OK', node)

        # Otherwise, recurse further down the tree
        violation, node.children[branch] = IntervalTree._insert_fix_and_update_max(
            node.children[branch], leaf)
        
        # We are unwinding the stack:

        # Check whether we have a violation. Skip step 2 if not.
        if violation == 'OK' or (violation == 'POTENTIAL' and not node.isRed):
            return ('OK', node)

        # Step 2: fix the violation where a red node has a red child.
        return IntervalTree._fix(node)

    @staticmethod
    def _fix(node):
        # Step 2: fix the violation where a red node has a red child.

        # If this node is red, we leave it to our parent to fix the violation.
        if node.isRed:
            return ('VIOLATION', node)

        # This node is black. Get all our red children to analyze the different
        # cases.
        redChildren = [ 
            (branch, child) for branch, child in enumerate(node.children)
            if IntervalTree._if_red(child) ]

        # We have to distinguish three cases:
        # Case 1 and 2: this node has one red child which has one red child.
        # Caes 3: this node has two red children (and one red grandchild).

        if len(redChildren) == 2:
            # We are in case 3. We change the color of this node and its two
            # children.
            for child in node.children:
                child.isRed = False
            node.isRed = True
            # This will fix the violation for the children, but might introduce
            # a violation where this node and its parent are both red.
            return ('POTENTIAL', node)

        # We are left with case 1 or 2. Get the one red child.
        branch, child = redChildren[0]

        # And look at one if its children to distinguish between case 1 and 2.
        # If our left child is red, get its right child.
        # If our right child is red, get its left child.
        grandChild = IntervalTree._if_red(child.children[1 - branch])
        if grandChild:
            # If that grand child is red, we are in case 2 which
            # requires an extra rotation to reduce to case 1.
            child.children[1 - branch] = grandChild.children[branch]
            child.update_max_value()
            grandChild.children[branch] = child
            node.children[branch] = grandChild
            child = grandChild
        # Case 1: perform rotation.
        node.children[branch] = child.children[1 - branch]
        node.update_max_value()
        node.isRed = True
        child.children[1 - branch] = node
        child.isRed = False

        # Case 1 and 2 have the violation resolved, no more fixing
        # required.
        return ('OK', child)

################################################################################
#
# TESTING

class _IntervalTreeTester(IntervalTree):
    """
    A test rig for IntervalTree. It will keep a separate plain list of all
    entries added by :py:meth:`insert`. It offers methods to compare the result
    of a brute force search in the plain list against an interval tree search.
    It also offers methods to do many consistency checks on the red black tree.

    Run the test rig::

        sage: _IntervalTreeTester.run_test()

    """

    def __init__(self):
        self._entries = []
        super(_IntervalTreeTester, self).__init__()

    def insert(self, interval, value):
        # Insert in plain list for reference
        self._entries.append((interval, value))
        # Insert into red black tree
        super(_IntervalTreeTester, self).insert(interval, value)

    def brute_force_find(self, interval):
        """
        Search plain list as ground truth to compare against.
        """
        return [ entry[1]
                 for entry in self._entries
                 if entry[0].overlaps(interval) ]

    def check_find_result(self, interval):
        if set(self.find(interval)) != set(self.brute_force_find(interval)):
            raise Exception("Different results: %r %r" % (
                    self.find(interval),
                    self.brute_force_find(interval)))

    def check_consistency(self):
        from sage.all import Infinity
        if self._root.isRed:
            raise Exception("Red root")
        _IntervalTreeTester._recursively_check_consistency(
            self._root, -Infinity, +Infinity)
        
    @staticmethod
    def _recursively_check_consistency(node, l, r):
        from sage.all import Infinity

        if not node:
            return -Infinity, 0

        if not (node.interval.lower() >= l):
            raise Exception("Node left  lower %r %r", node.interval.lower(), l)
        if not (node.interval.lower() <= r):
            raise Exception("Node right lower %r %r", node.interval.lower(), r)

        left_max, left_depth = (
            _IntervalTreeTester._recursively_check_consistency(
                node.children[LEFT], l, node.interval.lower()))

        right_max, right_depth = (
            _IntervalTreeTester._recursively_check_consistency(
                node.children[RIGHT], node.interval.lower(), r))
            
        if not max(left_max, right_max, node.interval.upper()) == node.max_value:
            raise Exception("Maximum incorrect")

        if left_depth != right_depth:
            raise Exception("Inconsistent black depths")

        if node.isRed:
            for child in node.children:
                if child and child.isRed:
                    raise Exception("Red node has red child")
        else:
            left_depth += 1
        
        return node.max_value, left_depth

    def print_tree(self):
        self.print_tree_recursively(self._root, 0)

    @staticmethod
    def print_tree_recursively(node, depth):
        if not node:
            return

        if not node.isRed:
            depth += 1

        _IntervalTreeTester.print_tree_recursively(node.children[0], depth)
            
        align = 6 * depth
        if node.isRed:
            align += 3

        print(align * " ", end = " ")
        if node.isRed:
            print("R", end = " ")
        else:
            print("B", end = " ")

        print(node.interval.lower())

        _IntervalTreeTester.print_tree_recursively(node.children[1], depth)
            
    @staticmethod
    def run_test():
        from sage.all import RIF, sin, Infinity
        
        intervals = [
            RIF(sin(1.2 * i), sin(1.2 * i) + sin(1.43 * i) ** 2)
            for i in range(200) ]

        t = _IntervalTreeTester()
        for i, interval in enumerate(intervals):
            t.insert(interval, i)
            if i % 50 == 0:
                for j in intervals:
                    t.check_find_result(j)
                t.check_consistency()

        num_true = len(intervals)
        num_have = len(t.find(RIF(-Infinity, Infinity)))
        if num_true != num_have:
            raise Exception("Inconsistent number of intervals: %d %d" % (
                    num_true, num_have))

def _doctest():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _doctest()
