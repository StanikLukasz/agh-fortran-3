# Fortran - project no. 3

Project in fortran for AGH UST classes.

### Zadanie

http://home.agh.edu.pl/~macwozni/fort/projekt3.pdf

### Teoria

* https://pl.wikipedia.org/wiki/Mno%C5%BCenie_macierzy
* https://pl.wikipedia.org/wiki/Metoda_eliminacji_Gaussa

### Praktyka

#### Specyfikacja

* U¿yty kompilator: ifort
* U¿yty jezyk: Fortran2008

#### Instrukcja uzytkowania

Aby skompilowac program nalezy uzyc odpowiedniej reguly make 
* `make main` - skompiluje program
* `make times` - uruchomi program 15 razy z rozmiarami macierzy od 100 do 1500 (co 100). Wynik zapisze w pliku tekstowym results.txt.

Po skompilowaniu mozna uruchomic program podajac rozmiar macierzy do przetworzenia w linii polecen. Program wywola cztery procedury:
* mnozenia macierzy sekwencyjnie
* mnozenia macierzy rownolegle
* eliminacji Gaussa sekwencyjnie
* eliminacji Gaussa rownolegle

Podczas wykonania program obliczy czasy wykonywania poszczegolnych procedur. Wyniki **dopisze** do pliku tekstowego results.txt.

W repozytorium znajduje sie rowniez prosta aplikacja napisana pythonie, ktora korzysta z procedur udostepnianych dzieki f2py przez sequential_lib.

#### Wyniki

Ponizej przedstawiam wykresy zaleznosci czasu obliczen od rozmiaru przetwarzanych macierzy  (w zakresie [100, 1500] o skali liniowej).

##### Mnozenie macierzy
###### tabela pomiarow
![table](https://github.com/StanikLukasz/agh-fortran-3/blob/master/times/mult_table.PNG)
###### wykres czasu od rozmiaru macierzy (czas obliczeñ podany w sekundach)
![chart](https://github.com/StanikLukasz/agh-fortran-3/blob/master/times/mult_chart.PNG)

##### Eliminacja Gaussa
###### tabela pomiarow
![chart](https://github.com/StanikLukasz/agh-fortran-3/blob/master/times/gauss_table.PNG)
###### wykres czasu od rozmiaru macierzy (czas obliczeñ podany w sekundach)
![chart](https://github.com/StanikLukasz/agh-fortran-3/blob/master/times/gauss_chart.PNG)

Wyraznie widac skrocenie czasu obliczen dzieki wykorzystaniu rownoleglosci. Efekt ten jest jednak malo widoczny (a czasem nawet odwrotny!) w przypadku danych wejsciowych o bardzo malych rozmiarach.