#Converstation Summary
Summary from Jan 9 2017 discussion between Dr Ben Evans and Xue Song regarding identification of splice variants.
- goal: identify total exons, max exons per transcript, min exons per transcript, location of variants(ends or middle)
- method: 
  - firstly, run Blast with set of assembled transcripts flagged with "no error" or "contain scaffold" against itself. 
  - Then used the Blast result and coordinate information from previous genomic blast to achieve goal
- critiria
  - transcripts that are just repeat region of each other does not count as splicing vairants. 
  - transcripts that identified as splicing variants should have genomic coordinate close to each other
 
