#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<complex>
using namespace std;
# define PI 3.14159265358979323846

// Evaluate the Transfer Function
complex<double> evaluateTransferFunction(const vector<double>& numeratorCoefficients, const vector<double>& denominatorCoefficients, complex<double> s)
{
    complex<double> numerator = 0.0;
    complex<double> denominator = 0.0;

    for(int i =0 ; i<numeratorCoefficients.size() ; ++i)
    {
        numerator += numeratorCoefficients[i] * pow(s, numeratorCoefficients.size()-i-1);
    }

    for(int i =0 ; i<denominatorCoefficients.size() ; ++i)
    {
        denominator += denominatorCoefficients[i] * pow(s, denominatorCoefficients.size()-i-1);
    }

    if(denominator == 0.0)
    {
        cout << "Error: Division by zero!" << endl;
        return 0.0;
    }
    return numerator/denominator;
}

// Calculate the Magnitude of the Transfer Function
double calculateMagnitude(const vector<double>& numeratorCoefficients, const vector<double> denominatorCoefficients, double frequency)
{
    complex<double> s = 1.0i * frequency;
    complex<double> transferFunction = evaluateTransferFunction(numeratorCoefficients, denominatorCoefficients, s);
    double magnitude = 20 * log10(abs(transferFunction));
    return magnitude;
}

// Calculate the Phase of the Transfer Function
double calculatePhase(const vector<double>& numeratorCoefficients, const vector<double> denominatorCoefficients, double frequency)
{
    complex<double> s = 1.0i * frequency;
    complex<double> transferFunction = evaluateTransferFunction(numeratorCoefficients, denominatorCoefficients, s);
    double phase = arg(transferFunction) * 180.0/PI;
    return phase;
}

// Start Gnuplot
void starteGnuplot(const string& dateiname)
{
    string kommando = "gnuplot " + dateiname;
    int ret = system(kommando.c_str());
    if(ret!=0)
    {
        cerr << "Fehler beim Ausfuhren von " << kommando;
    }
}

int main()
{
    double startFreq = 0.1;
    double endFreq = 10.0;
    int numPoints = 16;

    vector<double> freq(numPoints);
    vector<double> magnitude(numPoints);
    vector<double> phase(numPoints);
    
    ofstream magOutputFile("bode_data_magnitude.txt");
    ofstream phaseOutputFile("bode_data_phase.txt");
    const string plotAmplitude = "amplitude.gp"; 
    const string plotPhase = "phase.gp"; 

    double gainAtPhaseCrossover = 0.0;
    double phaseAtCrossover = 0.0;
    double frequencyAtPhaseCrossover = 0.0;
    double frequencyAtGainCrossover = 0.0;
    double gainMargin = 0.0;
    double phaseMargin = 0.0;

    double a,b,c,d,e,f,g,frequency;
    vector<double> numeratorCoefficients = {a,b,c};
    vector<double> denominatorCoefficients = {d,e,f,g};
    complex<double> s,result;
    vector<double> phaseProb(numPoints);
    vector<double> magnitudeProb(numPoints);


    int choice;
    do
    {
        cout << "MENU:\n";
        cout << "[1] Eingabe der Ubertragungsfunktion\n";
        cout << "[2] Plot der Amplitudengang\n";
        cout << "[3] Plot der Phasengang\n";
        cout << "[4] Berechnen der Amplitudenreserve und der Phasenreserve\n";
        cout << "[5] Beenden\n";
        cout << "Geben Sie die Zahl 1-5 ein\n";
        cin >> choice;

        switch (choice)
        {
            case 1: 
                // Enter the Coefficients for Transfer Function
                cout << "Enter coefficients for the numerator (a,b,c): ";
                cin >> a >> b >> c;
                cout << "Enter coefficients for the denominator (d,e,f,g): ";
                cin >> d >> e >> f >> g;

                numeratorCoefficients = {a,b,c};
                denominatorCoefficients = {d,e,f,g};

                // Print the Transfer Function
                cout << "H(s)= ";
                for(int i=0 ; i<numeratorCoefficients.size() ; ++i)
                {
                    cout << numeratorCoefficients[i] << "s^" << numeratorCoefficients.size()-i-1;
                    if(i< numeratorCoefficients.size()-1)
                        cout << " + ";
                } 
                cout << ") / (";
                for(int i=0 ; i<denominatorCoefficients.size() ; ++i)
                {
                    cout << denominatorCoefficients[i] << "s^" << denominatorCoefficients.size()-i-1;
                    if(i< denominatorCoefficients.size()-1)
                        cout << " + ";
                } 
                cout << ")" << endl;
                cout << "\n";

                // Only for Test Purpose. The Transfer Function will be evaluate at certain Point s.
                //cout << "Enter the value of s: ";
                //cin >> s;
                //result = evaluateTransferFunction(numeratorCoefficients, denominatorCoefficients, s);
                //cout << "\n";
                //cout << "Result for the expression for s= " << s << ": " << result << endl;
                break;
            case 2:
                // Plot Amplitudengang
                // Generate the Data Points
                for(int i=0 ; i<numPoints ; ++i)
                {
                    freq[i] = startFreq * pow(10, i * log10(endFreq/startFreq) / (numPoints - 1));
                    frequency = startFreq * pow(10, i * log10(endFreq/startFreq) / (numPoints-1));
                    magnitudeProb[i] = calculateMagnitude(numeratorCoefficients, denominatorCoefficients, frequency);
                }

                // Write into a File
                if(magOutputFile.is_open())
                {
                    for(int i=0 ; i<numPoints ; ++i)
                    {
                        magOutputFile   << freq[i] << "\t" << magnitudeProb[i]  << "\n" ;
                    }
                    magOutputFile.close();
                    cout << "\n" << "Data saved to 'bode_data_magnitude.txt'\n" << endl;
                }
                else
                {
                    cerr << "Unable to open the file for writing." << endl;
                    return 1;
                }

                // Plot the Amplitude with the Gnuplot
                starteGnuplot(plotAmplitude); 
                break;
            case 3:
                // Plot Phasengang 
                for( int i=0 ; i<numPoints ; ++i)
                {
                    freq[i] = startFreq * std::pow(10, i * log10(endFreq/startFreq) / (numPoints-1) );
                    frequency = startFreq * std::pow(10, i * log10(endFreq/startFreq) / (numPoints-1) );
                    phaseProb[i] = calculatePhase(numeratorCoefficients, denominatorCoefficients, frequency);
                }
                if (phaseOutputFile.is_open())
                {
                    for(int i=0 ; i<numPoints ; ++i)
                    {
                        phaseOutputFile << freq[i] << "\t" << phaseProb[i] << "\n" ;
                    }
                    phaseOutputFile.close();
                    cout << "\n" << "Data saved to 'bode_data_magnitude.txt'\n" << endl;
                }
                else
                {
                    std::cerr << "Unable to open the file for writing." << std::endl;
                    return 1;
                }
                starteGnuplot(plotPhase);
                break;
            case 4:
                // Amplitudenreserve und Phasenreserve
                for(int i=0 ; i<numPoints ; ++i)
                {
                    freq[i] = startFreq * std::pow(10, i * log10(endFreq/startFreq) / (numPoints-1) );
                    frequency = startFreq * std::pow(10, i * log10(endFreq/startFreq) / (numPoints-1) );
                    magnitudeProb[i] = calculateMagnitude(numeratorCoefficients, denominatorCoefficients, frequency);
                    phaseProb[i] = calculatePhase(numeratorCoefficients, denominatorCoefficients, frequency);

                    // Only for Test Purpose
                    //cout << freq[i] << "\t" << frequency << "\t" << magnitudeProb[i] << "\t" << phaseProb[i] <<endl;
                }
                
                // Calculate Gain Margin
                for(int i=0 ; i<numPoints ; ++i)
                {
                    if(phaseProb[i] <= -179.1)
                    {
                        gainAtPhaseCrossover = magnitudeProb[i];
                        frequencyAtPhaseCrossover = freq[i];
                        break;
                    }
                }

                // Calculate Phase Margin
                for(int i=0 ; i<numPoints ; ++i)
                {
                    if(magnitudeProb[i] <= 0.1)
                    {
                        phaseAtCrossover = phaseProb[i];
                        frequencyAtGainCrossover = freq[i];
                        break;
                    }
                }

                // Display Gain Margin
                gainMargin = - gainAtPhaseCrossover;
                cout << "Gain Margin: " << gainMargin << " dB\n";
                cout << "Frequency at Phase Crossover: " << frequencyAtPhaseCrossover << " rad/s" << endl;
                cout << "\n";
                // Display Phase Margin
                phaseMargin = phaseAtCrossover + 180.0;
                cout << "Phase Margin: " << phaseMargin << " degrees\n";
                cout << "Frequency at Gain Crossover: " << frequencyAtGainCrossover << " rad/s" << endl;
                cout << "\n";
                break;
            case 5:
                // Exit the Programm
                break;
            default:
                cout << "Falsche Eingabe. Probieren Sie noch einmal.\n" << endl;
        }

    } while (choice != 5);
            return 0;
    
}