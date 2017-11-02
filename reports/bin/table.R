################################################################################
# The effect of base modification on RNA polymerase and reverse transcriptase fidelity
# Copyright (C) 2017 New England Biolabs, Inc.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

args <- commandArgs(TRUE)

ifile  <- args[1]
prefix <- args[2]
type   <- args[3]

options(width=200)

df <- read.table(ifile, sep = ",", header = TRUE)

colsDesc <- c("SampleID","Enzyme","Template","Amplicon")
colsData <- c("AA","AC","AT","AG","CA","CC","CT","CG","TA","TC","TT","TG","GA","GC","GT","GG","Deletion","Insertion")

colsAT   <- c("AA","AC","AT","AG","TA","TC","TT","TG")
colsGC   <- c("CA","CC","CT","CG","GA","GC","GT","GG")

subsAT   <- c("AC","AT","AG","TA","TC","TG")
subsGC   <- c("CA","CT","CG","GA","GC","GT")

if ( type == "fse" ) {
    subsAT.rt <- c("rA:dG","rA:dA","rA:dC","rT:dT","rT:dG","rT:dC")
    subsGC.rt <- c("rC:dT","rC:dA","rC:dC","rG:dT","rG:dG","rG:dA")
} else {
    subsAT.rt <- c("dT:dC","dT:dT","dT:dG","dA:dA","dA:dC","dA:dG")
    subsGC.rt <- c("dG:dA","dG:dT","dG:dG","dC:dA","dC:dC","dC:dT")
}

## list of enzymes for analysis
enzymes <- c("ProtoScript II Reverse Transcriptase",
             "M-MuLV Reverse Transcriptase",
             "AMV Reverse Transcriptase",
             "Bst 2.0 DNA Polymerase",
             "Bst 3.0 DNA Polymerase")

## RNA templates (modification types)
templates <- c("Regular",
               "N6mA",
               "P-U",
               "5mC",
               "5mU",
               "5hmU")

## Amplicon names
amplicons <- c("DNA-1",
               "DNA-2",
               "DNA-3",
               "DNA-4")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Subset data before further processing
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## subset of samples with defined enzymes/templates
ix <- ( df$Enzyme %in% enzymes ) & ( df$Template %in% templates )

df0 <- df[ix,c(colsDesc,colsData,"D1","DX","I1","IX")]

## define factor levels for convenient sorting
df0$Enzyme   <- factor(df0$Enzyme,   levels = enzymes)
df0$Template <- factor(df0$Template, levels = templates)
df0$Amplicon <- factor(df0$Amplicon, levels = amplicons)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Aggregate by Enzyme,Template,Amplicon
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## combine counts for replicas
df1 <- aggregate(df0[c(colsData,"D1","DX","I1","IX")], df0[colsDesc], sum)
df1 <- df1[order(df1$Enzyme,df1$Template,df1$Amplicon),]

## derivative values
df1$Total                 <- apply(df1[c(colsAT,colsGC)],1,sum)

df1$ErrorSub              <- apply(df1[c(subsAT,subsGC)],1,sum)
df1$ErrorDelSingle        <- apply(df1[c("D1","DX")],1,sum)
df1$ErrorInsSingle        <- apply(df1[c("I1","IX")],1,sum)

df1$RateSub               <- df1$ErrorSub       / df1$Total
df1$RateDelSingle         <- df1$ErrorDelSingle / df1$Total
df1$RateInsSingle         <- df1$ErrorInsSingle / df1$Total

df1$RateSubSingleFraction <- df1$RateSub       / (df1$RateSub + df1$RateDelSingle + df1$RateInsSingle)
df1$RateDelSingleFraction <- df1$RateDelSingle / (df1$RateSub + df1$RateDelSingle + df1$RateInsSingle)
df1$RateInsSingleFraction <- df1$RateInsSingle / (df1$RateSub + df1$RateDelSingle + df1$RateInsSingle)

df1$RateTotal             <- df1$RateSub + df1$RateDelSingle + df1$RateInsSingle

write.table(df1,
            file = sprintf("%s-%s.1-all-samples.csv", prefix, type),
            append = FALSE,
            quote = TRUE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE )

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Aggregate by Enzyme,Template (collapse amplicons)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### compute average values and standard deviation
agg <- c("RateSub",
         "RateDelSingle",
         "RateInsSingle",
         "RateSubSingleFraction",
         "RateDelSingleFraction",
         "RateInsSingleFraction",
         "RateTotal"
         )

avg <- aggregate( df1[agg],        df1[c("Enzyme","Template")], mean   )
std <- aggregate( df1[agg],        df1[c("Enzyme","Template")], sd     )
tot <- aggregate( df1[c("Total")], df1[c("Enzyme","Template")], sum    )
num <- aggregate( df1[c("Total")], df1[c("Enzyme","Template")], length )
spe <- aggregate( df1[c(subsAT,subsGC,"ErrorSub")], df1[c("Enzyme","Template")], sum)

## rename column in std table to have unique column names
names(std) <- c("Enzyme","Template",gsub("Rate","Std",agg))
names(num) <- c("Enzyme","Template","Nsamples")
names(spe) <- c("Enzyme","Template",subsAT.rt,subsGC.rt,"ErrorSub")

## merge and sort
df2 <- merge( num, tot, by = c("Enzyme","Template") )
df2 <- merge( df2, avg, by = c("Enzyme","Template") )
df2 <- merge( df2, std, by = c("Enzyme","Template") )
df2 <- merge( df2, spe, by = c("Enzyme","Template") )

df2[,c(subsAT.rt,subsGC.rt)] <- df2[,c(subsAT.rt,subsGC.rt)] / df2$Total

df2 <- df2[order(df2$Enzyme,df2$Template),]

write.table(df2,
            file = sprintf("%s-%s.2-aggregate-by-template.csv", prefix, type),
            append = FALSE,
            quote = TRUE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE )

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Aggregate by Enzyme (collapse templates/amplicons)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

avg <- aggregate( df1[agg],        df1[c("Enzyme")], mean   )
std <- aggregate( df1[agg],        df1[c("Enzyme")], sd     )
tot <- aggregate( df1[c("Total")], df1[c("Enzyme")], sum    )
num <- aggregate( df1[c("Total")], df1[c("Enzyme")], length )
spe <- aggregate( df1[c(subsAT,subsGC,"ErrorSub")], df1[c("Enzyme")], sum)

## rename column in std table to have unique column names
names(std) <- c("Enzyme",gsub("Rate","Std",agg))
names(num) <- c("Enzyme","Nsamples")
names(spe) <- c("Enzyme",subsAT.rt,subsGC.rt,"ErrorSub")

## merge and sort
df3 <- merge( num, tot, by = c("Enzyme") )
df3 <- merge( df3, avg, by = c("Enzyme") )
df3 <- merge( df3, std, by = c("Enzyme") )
df3 <- merge( df3, spe, by = c("Enzyme") )

df3[,c(subsAT.rt,subsGC.rt)] <- df3[,c(subsAT.rt,subsGC.rt)] / df3$Total
## spe[,c(subsAT,subsGC)] <- spe[,c(subsAT,subsGC)] / spe[,"ErrorSub"]

df3 <- df3[order(df3$Enzyme),]

write.table(df3,
            file = sprintf("%s-%s.3-aggregate-by-enzyme.csv", prefix, type),
            append = FALSE,
            quote = TRUE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE )

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mutational spectrum
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df4 <- aggregate(df1[c(subsAT,subsGC,"ErrorSub","Total")], df1[c("Enzyme","Template")], sum)
df4 <- df4[order(df4$Enzyme,df4$Template),]

names(df4) <- c("Enzyme","Template",subsAT.rt,subsGC.rt,"ErrorSub","Total")

write.table(df4,
            file = sprintf("%s-%s.4-mutational-spectrum.csv", prefix, type),
            append = FALSE,
            quote = TRUE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE )

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Base specific rates
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df5 <- aggregate(df1[c(subsAT,subsGC,"ErrorSub","Total")],
                 df1[c("Enzyme","Template")],
                 sum)

df5["A"] <- apply(df5[c("AC","AG","AT")],1,sum)
df5["C"] <- apply(df5[c("CA","CG","CT")],1,sum)
df5["G"] <- apply(df5[c("GA","GC","GT")],1,sum)
df5["T"] <- apply(df5[c("TA","TC","TG")],1,sum)

df5[,c(subsAT,subsGC,"A","C","G","T")] <- df5[,c(subsAT,subsGC,"A","C","G","T")] / df5$Total

df5 <- df5[order(df5$Enzyme,df5$Template),]
    
names(df5) <- c("Enzyme","Template",subsAT.rt,subsGC.rt,"ErrorSub","Total","A","C","G","T")

write.table(df5,
            file = sprintf("%s-%s.5-base-specific-rates.csv", prefix, type),
            append = FALSE,
            quote = TRUE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE )
