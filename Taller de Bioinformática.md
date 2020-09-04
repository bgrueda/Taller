# Taller de Bioinformática

Tengo antecedentes de haber tomado algunos cursos de bioinformática (shell, perl, R, pyhton y otros básicos relacionados). Sin embargo con este nuevo proyecto que incluye un análisis de exoma aspiro a desarrollar mis habilidades en este campo y poder construir scripts más complejos por mí misma, de una forma más eficiente, así como resolver problemas ligados al desarrollo de los mismos. 


#### Algunos ejemplos de código que he realizado.  

**Perl**

```perl
#!/usr/bin/perl

# IDENTIFICA SI UNA SECUENCIA ES PALÍNDROME
print "\n Este programa te dice si tu secuencia es o no palíndrome \n";

print "\n ESCRIBE AQUÍ TU SECUENCIA: \n";
$rna = <STDIN>;
chomp ($rna);
$rna2 = reverse($rna);
if ($rna eq $rna2)
{print "Tu secuencia es palíndrome\n";}
else {print "Tu secuencia no es palíndrome\n";}

print "\n ¡Gracias por usar este programa! \n"
```

```perl
#!/usr/bin/perl

open (ABRIR, "Gene_sequence.txt");

while (<ABRIR>)
{
	@coord = split (/\t/);
	
		if (/(ECK[0-9]+)/)
		{
		$gen = $1; 
		}	
	%genes = ($gen, "$coord[2]-$coord[3]"); 
	foreach $k(sort keys%genes)
	{
	print "$k $genes{$k}\n";
	}
}

close(ABRIR);
```

```perl
#!/usr/bin/perl

use strict;
print "\n Este programa te dice qué tipo de triángulo tienes con base a las medidas de sus lados \n";

print "Cuanto mide el primer lado?\n";
my $a = <STDIN> ;

print "Cuanto mide el segundo lado?\n";
my $b = <STDIN> ;

print "Cuanto mide el tercer lado?\n";
my $c = <STDIN> ;

chomp ($a, $b, $c);

if ($a == $b and $b == $c)
{ print ("\n Ese es un triángulo EQUILÁTERO!!\n\n");}
elsif ($a == $b || $a == $c || $b == $c)
{ print ("\n  Ese es un triángulo ISÓSCELES!!\n\n");}
else 
{print ("\n Ese es un triángulo ESCALENO!!\n\n");}

```



**Pyhton**

```python
#!/usr/bin/python

import Tkinter

window = Tkinter.Tk()

res = open('tels.txt', 'a')


window.title('Widgets Example')

#window.wm_iconobotmap('alien.icp')

def callback():
	nam = nam.get()	
	print nam
	tel = tel.get()
	print tel
	re.write(nam + ', ' + tel + '\n')
	nam.delete(0,Tkinter,END)
	tel.delete(0,Tkinter,END)
	return 
	nam
	tel

def termina():
	res.close()
	window.quit()


lbl = Tkinter.Label(window, text='Datos', fg='purple')
lbl.pack()

Nam = Tkinter.Label(window, text= 'Nombre')
Nam.pack()
nam= Tkinter.Entry(window)
nam.pack()

Tel =Tkinter.Label(window, text='Telefono')
Tel.pack()
tel = Tkinter.Entry(window)
tel.pack()

proc = Tkinter.Button(window, text='PROCESA', command=callback)
proc.pack()
term = Tkinter.Button(window, text= 'SALIR', command=termina)
term.pack()

window.mainloop()
```



**R**

```R
setwd("~/Desktop/COVID/ROC")
library("pROC")

ROC <- read.table ("Ct.txt", 
                     header= T)
attach (ROC)

par(mfrow = c(2,4))

par(pty = "s")
roc(ROC$ID16, ROC$Ct16, plot=T, legacy.axes=T, percent=T
    , xlab="False Positive Percentage", ylab="True Postive Percentage",
    main= 'ROC 16.5 dilution', col="midnightblue", lwd=3, print.auc=T)

par(pty = "s")
roc(ROC$ID55, ROC$Ct55, plot=T, legacy.axes=T, percent=T
    , xlab="False Positive Percentage", ylab="True Postive Percentage",
    main= 'ROC 55 dilution', col="navy", lwd=3, print.auc=T)

par(pty = "s")
roc(ROC$ID82, ROC$Ct82, plot=T, legacy.axes=T, percent=T
    , xlab="False Positive Percentage", ylab="True Postive Percentage",
    main= 'ROC 82 dilution', col="dodgerblue4", lwd=3, print.auc=T)

par(pty = "s")
roc(ROC$ID165, ROC$Ct165, plot=T, legacy.axes=T, percent=T
    , xlab="False Positive Percentage", ylab="True Postive Percentage",
    main= 'ROC 165 dilution', col="royalblue3", lwd=3, print.auc=T)

par(pty = "s")
roc(ROC$ID330, ROC$Ct330, plot=T, legacy.axes=T, percent=T
    , xlab="False Positive Percentage", ylab="True Postive Percentage",
    main= 'ROC 330 dilution', col="blue3", lwd=3, print.auc=T)

par(pty = "s")
roc(ROC$ID550, ROC$Ct550, plot=T, legacy.axes=T, percent=T
    , xlab="False Positive Percentage", ylab="True Postive Percentage",
    main= 'ROC 550 dilution', col="blue2", lwd=3, print.auc=T)

par(pty = "s")
roc(ROC$ID1600, ROC$Ct1600, plot=T, legacy.axes=T, percent=T
    , xlab="False Positive Percentage", ylab="True Postive Percentage",
    main= 'ROC 1600 dilution', col="dodgerblue", lwd=3, print.auc=T)
```

```R
alfa = 0.05
de=3.6
n=100
me <- seq(14, 21, 0.1)
li <- 17.5 - qnorm(1-alfa/2)*(de/sqrt(100))
ls <- 17.5 + qnorm(1-alfa/2)*(de/sqrt(100))
b <- pnorm((ls-me)/(de/sqrt(n))) - pnorm((li-me)/(de/sqrt(n)))
pow <- 1-b

plot(me,pow, col='purple', lwd=2, pch=20, main= 'Gráfica de poder', xlab='Medias', ylab='Poder')
lines (me, pow, col='purple', lwd=2)
```





**Comandos básicos de shell y otras herramientas**

*1_EVALUACIÓN DE CALIDAD*

Cada archivo *.fastq se usa como input en FastQC

`fastqc`

Se remueven adaptadores y si es necesario se remueven nucleótidos de uno o ambos extremos de las lectuas

```shell
java -jar /route/trimmomatic-0.36.jar PE -phred33 file_R1.fastq file_R2.fastq file_R1_trimmed.fastq file_R1_unpaired.fastq file_R2_trimmed.fastq file_R2_unpaired.fastq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 CROP:150 HEADCROP:15 MINLEN:50
```

 Se vuelve a usar FastQC para asegurar que la calidad de las lecturas es adecuada.

Convertir de fastq a fasta

```shell
fastq_to_fasta -i file_R1_trimmed.fastq -o file_R1_trimmed.fa
```



*2_MAPEO* 

Mapear contra el genoma de referencia (X2, secuencia duplicada del genoma de referencia para abarcar las deleciones que crucen el origen) 

```shell
for file in *.fa; do ./blat chrM.fasta $file $file".psl"; done 
```

 

Concatenar archivo R1 y R2 para trabajar un solo archivo *.psl

```shell
tail -total_líneas-5 fileR2.psl > fileR2.psl 

cat fileR1.psl fileR2-5.psl > file.psl 
```



*3_IDENTIFICACIÓN DE DELECIONES*

Filtrar las deleciones 

```shell
for file in *.psl ; do sed 's/,/\t/g' $file | sed -E 's/\t\t/\t/g' | awk '($5 == 0 && $6 ==0 &&  $16 <=16569 && $18==2 && $19 >=15 && $20>=15)'| cut -f14-24 |sort -k1b,1 | uniq -c | sed 's/ *//' | sed 's/ /\t/' | cut -f2-12 | awk '$6+$7>=100 {print $0}' | awk '{print $10+$6-1, $11}' | sort -n > "breakpoints-"$file ; done 
```

 Corrección circular

```shell
awk '{if ($1>16569) print $1-16569, $2; else print $0}' breakpoints-file.psl | awk '{if ($2>16569) print $1, $2-16569; else print $0}' | sort -k 1n -k 2n > breakpoints-corregidos-file.txt
```

 Determinar frecuencias

```shell
breakpoints-corregidos-file.txt | sed 's/ *//' | sort -k1n -k2n > frecuencias-breakpoints.txt
```

Determinar tamaños de deleciones 

```shell
awk '{print $0, $3 - $2}' frecuencias-breakpoints.txt | awk '{if ($4 <= 0) print $1, $2, $3, $3+16569-$2 ; else print $0}' |  sort -k4n -k1n -k2n | sed 's/ /\t/g' > talla-frecuencias-breakpoints.txt
```



*4_SECUENCIAS REPETIDAS DIRECTAS*

Generación de breakpoints aleatorios para determinar la probabilidad de asociación con secuencias repetidas directas

​	**En R**

```R
x = sample (1:16569, 20000, replace=T)

y = sample (1:16569, 20000, replace=T)

tabla <- data.frame (x,y)

write.csv (tabla, "DR.csv") 
```

 Añadir al archivo de breakpoints la columna “chrM” y un identificador a cada evento: **RD-file.csv**

Crear un archivo para las coordenadas de start y otro para las coordenadas de end

```shell
awk '{print $1"\t" $2-20"\t" $2"\t" $5"\t" $4"\t"}' RD-file.csv | grep -v "-" > RD-file-start.csv 

awk '{print $1"\t" $3-20"\t" $3"\t" $5"\t" $4"\t"}' RD-file.csv | grep -v "-" > RD-file-end.csv
```

Usar Bedtools para recuperar las secuencias con base a las coordenadas

```shell
bedtools getfasta -name -tab -fi chrM.fasta -bed RD-file-start.csv -fo RD-start-intermediate.csv ; sed 's/:/ /g' RD-start-intermediate.csv > RD-file-start-sequences.csv

bedtools getfasta -name -tab -fi chrM.fasta -bed RD-file-end.csv -fo RD-end-intermediate.csv ; sed 's/:/ /g' RD-end-intermediate.csv >  RD-file-end-sequences.csv
```

 Alinear de acuerdo a los identificadores (Archivos separados por \t)

```shell
join RD-file-start.csv RD-file-end.csv > join-file1.txt

awk '{print$4"\t" $2"\t" $3}' RD-file.csv > join-file2.txt

join join-file1.txt join-file2.txt | awk '{print $1"\t" $4"\t" $5"\t" $2"\t" $3}' > RD-file-ready.csv
```

En este momento es posible proceder a graficar. 
