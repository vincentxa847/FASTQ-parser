import argparse

parser = argparse.ArgumentParser()
parser.add_argument('fastq', nargs='*') # input more than one files
parser.add_argument('--quality_filter', required=False, type=int, help='for_average_quality')
args = parser.parse_args()


class FASTQ:
    def __init__(self, name):
        self.information = []  # store sequence object
        self.name = name
        self.linenumber = 0

    def __str__(self):
        return f'FASTQ name = {self.name} '

    def Add_sequence(self):  # add the sequence object into FASTQ file object
        self.information.append(SequenceObject)

    def storeFastqtoList(self):
        lineWithoutNewlineCharacter = []

        with open(self.name, 'r') as fastq:  # open to store all the lines in a list and count the number of lines
            lineInFile = fastq.readlines()

        self.linenumber = (len(lineInFile))
        for i in lineInFile:
            lineWithoutNewlineCharacter.append(i.rstrip()) # remove the newline character
        return lineWithoutNewlineCharacter

class BLOCK:
    def __init__(self, i, currentBlock):
        self.file = i
        self.block = currentBlock
        self.id = ''
        self.sequence = ''
        self.ThirdLine = ''
        self.quality = ''
        self.bad = False
        self.GC_Content_percent = 0
        self.NumberOfMissingBase = 0
        self.total_score = 0

    def __str__(self):
        return f'sequence in {self.file}'

    def block(self):
        return self.block

    def resetBlockCondition(self):
        self.bad = False

    def blockInformationPrintAndResetValue(self):
        if self.bad:
            pass
        else:
            print(f"{self.file} {self.id} {self.NumberOfMissingBase} {self.GC_Content_percent}")
        self.NumberOfMissingBase = 0
        self.GC_Content_percent = 0
        self.total_score = 0


    def fastawrite(self):
        if (self.bad == False) and (self.total_score > int(args.quality_filter)):
            FASTA.write(f'> {self.id}\n{self.sequence}\n')


    def identifier(self):
        Firstline = lineInFile[self.block]
        identifier = []
        for k in Firstline:  # take the identifier from line 1
            identifier.append(k)
            if k == " ":
                break

        if Firstline == '':
            # an empty line here, still need to count because it not jumps into next line
            # the condition for empty line must be in the first if statement because empty line cannot be indexed ( will cause index error )
            self.bad = True
            self.block = self.block + 1
            print("{} no id in this entry".format(self.file))
            # raise Exception("mark1")

        elif ("@" != identifier[0]) and (characterForLine2  not in Firstline) and ("+" not in identifier) and (characterForLine4 not in identifier): # for bad entry
            # not a typical id line format, so not a good entry here
            # need to make sure "@" != identifier[0] not because the reading frame shift to other line ( id line loss )
            self.bad = True
            self.block = self.block + 1
            print("{}  id in bad format".format(self.file))
            # raise Exception("mark2") # using raise exception to check error

        elif "@" == identifier[0]: # normal line1
            identifier.pop(0)
            self.id = ''.join(identifier)  # store the FASTQ entry line 1
            self.block = self.block + 1
            # raise Exception("mark3")

        elif (identifier.count("A") >= 5) or (identifier.count("T") >= 5) or (identifier.count("C") >= 5) or (identifier.count("G") >= 5) or ("+" in identifier) or (characterForLine4 in identifier):
            # in a normal situation A,T,C,G would not appear in id for more than 5 times, so I use this as a condition to distinguish identifier from sequence
            # also "+" and characterforline4 should not appear in id too
            self.bad = True # line 1 disappear and directly go to other line
            print("{} no id in this entry".format(self.file))
            #raise Exception("mark4")


    def missing_base_AND_GC_content(self):
        Secondline = lineInFile[self.block]
        listForSecondline = [] # create a list here is for distinguishing line 4 (character for line4 in list)
        for v in Secondline:
            listForSecondline.append(v)
        if Secondline == '':
        # an empty line here, still need to count because it not jumps into next line
            if not self.bad:
                self.bad = True
                self.block = self.block + 1
                print("{} {} no line2 in this entry".format(self.file, self.id))
                # raise Exception("mark7")
            else:
                self.block = self.block + 1
                self.bad = True
                # raise Exception("mark8")

        elif (characterForLine2 not in Secondline) and ("@" != Secondline[0]) and ("+" != Secondline[0]) and (characterForLine4 not in listForSecondline):
            # for bad entry
            # not a typical sequence line, so not a good entry here
            # need to make sure no A,T,C,G, here is not because the reading frame shift to other line ( sequence line loss )
            if not self.bad:
                print("{} {} line2 in bad format or in lower case".format(self.file, self.id))
                self.bad = True
                self.block = self.block + 1
                #raise Exception("mark5")
            else:
                self.block = self.block + 1
                self.bad = True
                #raise Exception("mark6")

        elif characterForLine2 in Secondline:
            # normal line2 here
            if self.bad:
                self.block = self.block + 1
                #raise Exception("mark9")
                # if the id line already bad, don't print just count
            else:
                self.block = self.block + 1
                self.sequence = Secondline  # store the FASTQ entry line 2
                NumberOfMissingBase = 0
                GC_Content = 0
                for i in Secondline:
                    if i == "N":
                        NumberOfMissingBase += 1
                    if i == "G":
                        GC_Content += 1
                    if i == "C":
                        GC_Content += 1
                GC_Content_percent = (GC_Content / len(Secondline)) * 100
                self.GC_Content_percent += GC_Content_percent
                self.NumberOfMissingBase += NumberOfMissingBase
                #raise Exception("mark10")

        elif "@" == Secondline[0] or "+" == Secondline[0] or characterForLine4 in listForSecondline:
            # line2 loss and go down to other line
            if not self.bad:
                print("{} {} no line2 in this entry".format(self.file,self.id))
                self.block = self.block + 0 # line 2 total lost and go down to next line
                self.bad = True
                # raise Exception("mark11")
            else:
                self.bad = True
                #raise Exception("mark12")

    def line3(self):
        Thirdline = lineInFile[self.block]  # skip the third line (no important information here)
        listForline3 = []
        for l in Thirdline:
            listForline3.append(l)
        if Thirdline == '':
            # an empty line here, still need to count because it not jumps into next line
            if not self.bad:
                self.bad = True
                self.block = self.block + 1
                print("{} {} no line3 in this entry".format(self.file,self.id))
                #raise Exception("mark13")
            else:
                self.block = self.block + 1
                self.bad = True
                #raise Exception("mark14")

        elif ('+' != Thirdline[0]) and (characterForLine2 not in Thirdline) and ("@" != Thirdline[0]) and (characterForLine4 not in listForline3):  # for bad entry
            # not a typical thirdline,  so not a good entry here
            # need to make sure that not start with + is not because the reading frame shift to other line (thirdline loss)
            if not self.bad:
                self.bad = True
                print("{} {} line 3 in bad format".format(self.file, self.id))
                self.block = self.block + 1
                #raise Exception("mark15")
            else:
                self.block = self.block + 1
                #raise Exception("mark16")

        elif '+' == Thirdline[0]:
            # normal entry
            self.block = self.block + 1
            self.ThirdLine = Thirdline  # store the FASTQ entry line 3
            # raise Exception("mark17")

        elif "@" or characterForLine2 or characterForLine4 in Thirdline:
            # skip to other line
            if not self.bad:
                print("{} {} no line3 in this entry".format(self.file,self.id))
                self.bad = True
                #raise Exception("mark18")
            else:
                self.bad=True
                #raise Exception("mark19")

    def quality(self):
        Forthline = lineInFile[self.block]
        if Forthline == '':
            # an empty line here, still need to count because it not jumps into next line
            if not self.bad:  # Already a bad entry, so just count
                self.bad = True
                self.block = self.block + 1
                print("{} {}no line4 in this entry".format(self.file, self.id))
                # raise Exception("mark20")
            else:
                self.block = self.block + 1
                self.bad = True
                # raise Exception("mark21")

        for character in characterForLine4:
            if character in Forthline:
                #raise Exception("mark22")
                self.block = self.block + 1
                self.quality = Forthline  # store the FASTQ entry line 4
                total_score = 0
                for i in Forthline:  # take the quality value from forth line
                    score = ord(i)
                    total_score += score
                self.total_score += total_score
                break
        else:
            if not self.bad:
                print("{} {}no line4 in this entry".format(self.file, self.id))
                self.bad = True
            else:
                self.bad = True

characterForLine2 = "A" and "T" and "C" and "G" and "AA" and "TT" and "CC" and "GG"
characterForLine2LowerCase = characterForLine2.lower()
# a normal sequence will have equal number of A,T,C,G, just in case the entry sequence come from special region, so put additional condition here
characterForLine4 = ["%", "*" , "?" , "!" , "(" , ")" , "&" , "<" , ">" , "," ]
# common symbols in quality score, which should only appear in line4


for i in args.fastq:
    if ".fastq" not in i:  # check whether the input file is a fastq file
        print("{} not a fastq file".format(i))
        continue
    try:
        FileObject = FASTQ(i)  # for each fastaq file, create a file object
        fastq = open(i, 'r')  # reopen to store all the lines in a list
        lineInFile = FASTQ.storeFastqtoList(self=FileObject) # count line number here
        FASTAName = str(i).replace('fastq', 'fasta')
        FASTA = open(FASTAName, 'w')
        currentBlock = 0
        SequenceObject = BLOCK(i, currentBlock)
        try:
            while currentBlock <= FileObject.linenumber:  # iterate through every line in file
                BLOCK.identifier(self=SequenceObject)
                BLOCK.missing_base_AND_GC_content(self=SequenceObject)
                BLOCK.line3(self=SequenceObject)
                BLOCK.quality(self=SequenceObject)
                FASTQ.Add_sequence(self=FileObject)
                currentBlock = int(BLOCK.block(self=SequenceObject))
                BLOCK.fastawrite(self=SequenceObject)
                BLOCK.blockInformationPrintAndResetValue(self=SequenceObject)
                BLOCK.resetBlockCondition(self=SequenceObject)
                if currentBlock >= FileObject.linenumber:
                    break
            FASTA.close()
        except IndexError:
            print(f"{i} may lead to infinite loop")
    except FileNotFoundError:
        print(f"file {i} not found")