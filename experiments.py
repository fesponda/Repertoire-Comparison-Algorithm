
# takes two files, generates n samples of size m from each

num_experiments = 20
def generate_samples(df1,df2,id=''):
    l1=len(df1)
    l2=len(df2)
    min=l1
    if l2<min:
        min=l2
    sampleSize=int(0.8*min)

    arch = 'Spleen_1'
    for i in range(num_experiments):
        name = arch + 'exp_' + str(i) + id

        df = df1.sample(sampleSize)
        df = df.sort_values(by=['count (templates/reads)'], ascending=False)
        df.to_csv(dataPath + name + 'Sample.tsv', sep='\t', index=False)

        df = df2.sample(sampleSize)
        df = df.sort_values(by=['count (templates/reads)'], ascending=False)
        df.to_csv(dataPath + name + 'Sample.tsv', sep='\t', index=False)



