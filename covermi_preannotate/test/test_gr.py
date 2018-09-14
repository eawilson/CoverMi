import unittest
from gr import Gr, Entry, Chrom, Iterate, bisect_left



class TestEntry(unittest.TestCase):

    def test_constructor(self):
        self.assertEqual(Entry(1, 10, 20), Entry("chr1", 10, 20))
        self.assertEqual(Entry(23, 10, 20), Entry("23", 10, 20))
        self.assertEqual(Entry(23, 10, 20), Entry(23, 10, 20))
        self.assertEqual(Entry(23, 10, 20), Entry("chrX", 10, 20))
        self.assertEqual(Entry(23, 10, 20), Entry("X", 10, 20))
        self.assertEqual(Entry(23, 10, 20), Entry("X", 10, 20, "?"))
        self.assertEqual(Entry(23, 10, 20), Entry("X", 10, 20, "?", "+"))
        self.assertEqual(Entry(23, 10, 20).name, ".")
        self.assertEqual(Entry(23, 10, 20).strand, ".")
        with self.assertRaises(KeyError):
            Entry(0, 10, 20)
            Entry(26, 10, 20)
            Entry("?", 10, 20)

    def test_repr(self):
        self.assertEqual(eval(repr(Entry(23, 10, 20))), Entry(23, 10, 20,))
        self.assertEqual(repr(Entry(23, 10, 20, "?", "+")), 'Entry(Chrom(23), 10, 20, "?", "+")')

    def test_location(self):
        self.assertEqual(Entry(23, 10, 20).location, "chrX:10-20")

    def test_eq_hash(self):
        self.assertNotEqual(Entry(23, 10, 20), set([Entry(23, 10, 21)]))
        self.assertNotEqual(Entry(23, 10, 20), set([Entry(23, 9, 2)]))
        self.assertNotEqual(Entry(23, 10, 20), set([Entry(21, 10, 21)]))
        self.assertIn(Entry(23, 10, 20, "?", "+"), set([Entry(23, 10, 20)]))
        self.assertNotIn(Entry(23, 10, 20, "?", "+"), set([Entry(23, 11, 20, "?", "+")]))


class TestIterate(unittest.TestCase):

    def test_bisect_left(self):
        a = Entry(1, 5, 9)
        b = Entry(1, 10, 14)
        c = Entry(1, 20, 24)
        d = Entry(1, 30, 34)
        e = Entry(1, 40, 44)
        array = Gr([e, d, c, b, a]).data[1]
        self.assertEqual(bisect_left(array, 30, 0, len(array)), 3)
        array = Gr([e, d, c, b, a, Entry(1, 100, 130)]).data[1]
        self.assertEqual(bisect_left(array, 30, 0, len(array)), 0)
        array = Gr([e, d, c, b, a, Entry(1, 100, 120)]).data[1]
        self.assertEqual(bisect_left(array, 30, 0, len(array)), 1)
        array = Gr([e, d, c, b, a, Entry(1, 100, 200)]).data[1]
        self.assertEqual(bisect_left(array, 30, 0, len(array)), 0)
        array = Gr([e, d, c, b, a, Entry(1, 100, 120)]).data[1]
        self.assertEqual(bisect_left(array, 400, 0, len(array)), len(array))

    def test_iterate(self):
        a = Entry(1, 5, 9)
        b = Entry(1, 10, 14)
        c = Entry(1, 20, 24)
        d = Entry(1, 30, 34)
        e = Entry(1, 40, 44)
        array = Gr([e, d, c, b, a]).data[1]
        self.assertEqual(list(Iterate(array).yield_overlapping(1, 4)), [])
        self.assertEqual(list(Iterate(array).yield_overlapping(1, 5)), [a])
        self.assertEqual(list(Iterate(array).yield_overlapping(24, 30)), [c, d])
        self.assertEqual(list(Iterate(array).yield_overlapping(24, 40)), [c, d, e])
        self.assertEqual(list(Iterate(array).yield_overlapping(24, 1000)), [c, d, e])
        array = Gr([e, d, c, b, a, Entry(1, 1000, 2000)]).data[1]
        self.assertEqual(list(Iterate(array).yield_overlapping(29, 35)), [d])
        x = Entry(1, 9, 30)
        array = Gr([e, d, c, b, a, x]).data[1]
        self.assertEqual(list(Iterate(array).yield_overlapping(30, 30)), [x, d])


class TestGr(unittest.TestCase):

    def test_repr(self):
        gr = Gr([Entry(23, 10, 20), Entry(1, 56, 88, "?", "+")])
        self.assertEqual(eval(repr(gr)), gr)
        self.assertEqual(repr(gr), 'Gr([Entry(Chrom(1), 56, 88, "?", "+"), Entry(Chrom(23), 10, 20, ".", ".")])')

    def test_initialisation(self):
        a = Entry(1, 5, 9, "c")
        b = Entry(1, 10, 20, "d")
        c = Entry(1, 10, 21, "c")
        d = Entry(1, 30, 40, "a")
        e = Entry(2, 10, 20, "b")
        gr = Gr([e, d, c, b, a])
        self.assertEqual(gr, Gr([a, b, c, d, e]))
        self.assertEqual(list(gr), [a, b, c, d, e])

    def test_properties(self):
        a = Entry(1, 5, 9, "c")
        b = Entry(1, 10, 20, "d")
        c = Entry(1, 10, 21, "c")
        d = Entry(1, 30, 40, "a")
        e = Entry(2, 10, 20, "b")
        gr = Gr([e, d, c, b, a])
        self.assertEqual(gr.number_of_components, 5)
        self.assertEqual(gr.bases, 50)
        self.assertEqual(gr.is_empty, False)
        self.assertEqual(Gr().is_empty, True)
        self.assertEqual(gr.names, ["a", "b", "c", "d"])
        self.assertEqual(gr.names_as_string, "a, b, c, d") # Need to test exons
        self.assertEqual(gr.locations_as_string, "chr1:5-9, chr1:10-20, chr1:10-21, chr1:30-40, chr2:10-20")
        self.assertEqual(list(gr.sorted_by_name), [d, e, a, c, b])
                     
    def test_hash_eq(self): # test with variants
        a = Entry(1, 5, 9, "c")
        b = Entry(1, 10, 20, "d")
        c = Entry(1, 10, 21, "c")
        d = Entry(1, 30, 40, "a")
        e = Entry(2, 10, 20, "b")
        gr = Gr([e, d, c, b, a])
        self.assertEqual(gr, Gr([a, b, c, d, e]))
        self.assertNotEqual(gr, Gr([b, c, d, e]))
        self.assertNotEqual(gr, Gr([a, a, b, c, d, e]))
        self.assertIn(gr, set([Gr([a, b, c, d, e])]))
        self.assertNotIn(gr, set([Gr([b, c, d, e])]))
        self.assertNotIn(gr, set([Gr([a, a, b, c, d, e])]))

    def merged(self):
        a = Entry(1, 5, 9, "c")
        b = Entry(1, 10, 20, "d")
        c = Entry(1, 11, 15, "c")
        d = Entry(1, 17, 40, "a")
        e = Entry(1, 42, 60, "b")
        f = Entry(2, 42, 60, "e")
        self.assertEqual(Gr([a, a, a]).merged, Gr([a]))
        self.assertEqual(Gr([a, b, c, d, e, f]).merged, Gr([Entry(1, 5, 40), e, f]))
        self.assertEqual(Gr([b, e]).merged, Gr([b, e]))
        self.assertEqual(Gr([b, c]).merged, Gr([b]))

    def test_inverted(self):
        a = Entry(1, 5, 9, "c")
        b = Entry(1, 10, 20, "d")
        c = Entry(1, 10, 21, "c")
        d = Entry(1, 30, 40, "a")
        e = Entry(2, 10, 20, "b")
        gr = Gr([a])
        self.assertEqual(gr.inverted.inverted, gr)
        self.assertEqual(gr.inverted.number_of_components, 26)
        self.assertEqual(len(set(entry.chrom for entry in gr.inverted)), 25)
        gr = Gr([e, d, c, b, a])
        self.assertEqual(gr.inverted.inverted, gr.merged)
        self.assertEqual(gr.inverted.number_of_components, 28)

    def test_overlapped_by(self):
        a = Entry(1, 5, 9)
        b = Entry(1, 10, 20)
        c = Entry(1, 11, 15)
        d = Entry(1, 17, 40)
        e = Entry(1, 42, 60)
        f = Entry(2, 42, 60)
        gr = Gr([a, b, c, d, e, f])
        gr2 = Gr([a, b, c, d, e, f]*2)
        x = Entry(1, 12, 13)
        self.assertEqual(gr.overlapped_by(Gr([Entry(1,2,3)])), Gr())
        self.assertEqual(gr.overlapped_by(gr), gr)
        self.assertEqual(gr.overlapped_by(gr2), gr)
        self.assertEqual(gr2.overlapped_by(gr), gr2)
        self.assertEqual(gr2.overlapped_by(gr2), gr2)
        self.assertEqual(gr.overlapped_by(Gr([x])), Gr([x]*2))
        self.assertEqual(gr.overlapped_by(Gr([Entry(1, 9, 10)])), Gr([Entry(1, 9, 9), Entry(1, 10, 10)]))
        self.assertEqual(Gr([Entry(1, 2, 200)]).overlapped_by(gr), Gr([Entry(1, 5, 40), e]))

    def test_touched_by(self):
        a = Entry(1, 5, 9)
        b = Entry(1, 10, 20)
        c = Entry(1, 11, 15)
        d = Entry(1, 17, 40)
        e = Entry(1, 42, 60)
        f = Entry(2, 42, 60)
        gr = Gr([a, b, c, d, e, f])
        self.assertEqual(gr.touched_by(Gr([Entry(1, 3, 4)])), Gr())
        self.assertEqual(gr.touched_by(Gr([Entry(1, 5, 5)])), Gr([a]))
        self.assertEqual(gr.touched_by(Gr([Entry(1, 1, 100)])), Gr([a, b, c, d, e]))
        self.assertEqual(gr.touched_by(Gr([Entry(1, 16, 17)])), Gr([b, d]))
        self.assertEqual(Gr([b]+[c]*1000+[d]).touched_by(Gr([Entry(1, 19, 19)])), Gr([b, d]))

        self.assertEqual(gr.touched_by(Gr([Entry(1, 3, 4)]*2)), Gr())
        self.assertEqual(gr.touched_by(Gr([Entry(1, 5, 5)]*2)), Gr([a]))
        self.assertEqual(gr.touched_by(Gr([Entry(1, 1, 100)]*2)), Gr([a, b, c, d, e]))
        self.assertEqual(gr.touched_by(Gr([Entry(1, 16, 17)]*2)), Gr([b, d]))
        self.assertEqual(Gr([b]+[c]*1000+[d]).touched_by(Gr([Entry(1, 19, 19)]*2)), Gr([b, d]))

    def test_not_touched_by(self):
        a = Entry(1, 5, 9)
        b = Entry(1, 10, 20)
        c = Entry(1, 11, 15)
        d = Entry(1, 17, 40)
        e = Entry(1, 42, 60)
        f = Entry(2, 42, 60)
        gr = Gr([a, b, c, d, e, f])
        self.assertEqual(gr.not_touched_by(Gr([Entry(1, 3, 4)])), gr)
        self.assertEqual(gr.not_touched_by(Gr([Entry(1, 5, 5)])), Gr([b, c, d, e, f]))
        self.assertEqual(gr.not_touched_by(Gr([Entry(1, 1, 100)])), Gr([f]))
        self.assertEqual(gr.not_touched_by(Gr([Entry(1, 16, 17)])), Gr([a, c, e, f]))
        self.assertEqual(Gr([b]+[c]*1000+[d]).not_touched_by(Gr([Entry(1, 19, 19)])), Gr([c]*1000))

        self.assertEqual(gr.not_touched_by(Gr([Entry(1, 3, 4)]*2)), gr)
        self.assertEqual(gr.not_touched_by(Gr([Entry(1, 5, 5)]*2)), Gr([b, c, d, e, f]))
        self.assertEqual(gr.not_touched_by(Gr([Entry(1, 1, 100)]*2)), Gr([f]))
        self.assertEqual(gr.not_touched_by(Gr([Entry(1, 16, 17)]*2)), Gr([a, c, e, f]))
        self.assertEqual(Gr([b]+[c]*1000+[d]).not_touched_by(Gr([Entry(1, 19, 19)]*2)), Gr([c]*1000))

    def test_subranges_covered_by(self):
        a = Entry(1, 5, 9)
        b = Entry(1, 10, 20)
        c = Entry(1, 11, 15)
        d = Entry(1, 17, 40)
        e = Entry(1, 42, 60)
        f = Entry(2, 42, 60)
        x = Entry(1, 16, 30)
        gr = Gr([a, x])
        self.assertEqual(gr.subranges_covered_by(Gr([b])), Gr())
        self.assertEqual(gr.subranges_covered_by(Gr([b, c, d])), Gr([x]))
        self.assertEqual(gr.subranges_covered_by(Gr([b]*2)), Gr())
        self.assertEqual(gr.subranges_covered_by(Gr([b, c, d]*2)), Gr([x]))

    def test_combined_with(self):
        a = Entry(1, 5, 9)
        b = Entry(1, 10, 20)
        c = Entry(1, 11, 15)
        d = Entry(1, 17, 40)
        e = Entry(1, 42, 60)
        gr = Gr([a, b, c, d, e])
        self.assertEqual(gr.combined_with(gr), Gr([a, b, c, d, e]*2))
        self.assertEqual(Gr([a]).combined_with(Gr([b])), Gr([a, b]))

    def test_subset(self): # need to test exons
        a = Entry(1, 5, 9, "x")
        b = Entry(1, 10, 20, "b")
        c = Entry(1, 10, 21, "x")
        d = Entry(1, 30, 40, "d")
        e = Entry(2, 10, 20, "x")
        gr = Gr([a, b, c, d, e])
        self.assertEqual(gr.subset(["x"]), Gr([a, c, e]))





if __name__ == '__main__':
    unittest.main()
