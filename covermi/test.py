from covermi import Gr, Entry


one = Gr([Entry("chr1", 10, 20),
          Entry("chr1", 30, 40),
          Entry("chr1", 35, 45),
          Entry("chr2", 10, 20),
          Entry("chr2", 30, 40),
          Entry("chr2", 35, 45)])

print(one.touched_by(Entry("chr1", "5", "9")))


