library(Seurat)
#https://www.jianshu.com/p/5b26d7bc37b7


getwd()
path="G:/silicosis/geo/GSE184854_scRNA-seq_mouse_lung_ccr2/GSE184854_RAW/"
setwd(path)
getwd()


#new_counts=read.table(file = "G:/silicosis/geo/GSE184854_scRNA-seq_mouse_lung_ccr2/GSE184854_RAW/GSM5598265_matrix_inflection_demulti_SilicaWT.txt/GSM5598265_matrix_inflection_demulti_SilicaWT.txt")
head(new_counts)
#鎴栬€呭彲浠ョ洿鎺? create鑹瞮rat object锛坈ounts=counts锛?
mydata <- CreateSeuratObject(counts = read.table(file = "G:/silicosis/geo/GSE184854_scRNA-seq_mouse_lung_ccr2/GSE184854_RAW/GSM5598265_matrix_inflection_demulti_SilicaWT.txt/GSM5598265_matrix_inflection_demulti_SilicaWT.txt"),
                             min.cells = 10, project = "mydata_scRNAseq")

dim(mydata)
Idents(mydata)
mydata$orig_stim=Idents(mydata)
table(Idents(mydata))
colnames(mydata)
library(stringr)

teststring=colnames(mydata)[1:50]
str_split(teststring,"_",simplify = T)
str_split(teststring,"_",simplify = T)[,2]

table(str_split(colnames(mydata),"_",simplify = T)[,2])
table(str_split(colnames(mydata),"_",simplify = T)[,1])


mydata$percent.mt=PercentageFeatureSet(mydata,pattern = "^Mt")
Idents(mydata)


grepl(colnames(mydata),pattern = 'A0301')
colnames(mydata)[grepl(colnames(mydata),pattern = 'A0301')]

table(Idents(mydata))
mydata=RenameIdents(mydata,'A0301'='day_21_1','A0302'='day_21_2','A0303'='day_21_3',
                    'A0304'='day_7_1','A0305'='day_7_2','A0306'='day_7_3',
                    'A0307'='day_3_1','A0308'='day_3_2','A0309'='day_3_3',
                    'doublet'='doublet')
mydata$mystim=Idents(mydata)
table(mydata$mystim)
dim(mydata)


Idents(mydata)=mydata$orig_stim
mydata=RenameIdents(mydata,'A0301'='Day_21','A0302'='Day_21','A0303'='Day_21',
                    'A0304'='Day_7','A0305'='Day_7','A0306'='Day_7',
                    'A0307'='Day_3','A0308'='Day_3','A0309'='Day_3',
                    'doublet'='Doublet')
mydata$my_stim_3=Idents(mydata)


save(mydata,file = "G:/silicosis/geo/GSE184854_scRNA-seq_mouse_lung_ccr2/GSE184854_RAW/mydata.rds")
