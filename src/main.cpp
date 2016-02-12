#include <include/stdafx.h>
#include <include/mrvm.h>
#include <QtCore/QCoreApplication>
#include <iostream>
#include <cstdio>
#include <QFileInfo>
#include <arrayfire.h>

using namespace std;
using namespace af;


int main(int argc, char *argv[])
{
    bool hasCheckedArgs = false;

    while (true)
    {
        std::vector<QString> arguments;

        if(!hasCheckedArgs && argc > 1)
        {
            for(int i = 1; i < argc ; i++)
            {
                arguments.push_back(QString(argv[i]));
            }

            hasCheckedArgs = true;
        }
        else
        {
            cout << "Enter command: " ;
            std::string line;
            std::getline(std::cin, line);

            cout << endl;

            std::istringstream liness(line);

            std::string command;

            while (liness >> command)
            {
                QString qcommand = QString::fromStdString(command);

                if(qcommand[0] == '\"' || qcommand[0] == '\'')
                {
                    qcommand = qcommand.remove(0,1);
                }

                if(qcommand[qcommand.length() -1] == '\"' || qcommand[qcommand.length() -1] == '\'')
                {
                    qcommand = qcommand.remove(qcommand.length()-1,1);
                }

                arguments.push_back(qcommand);
            }

            cout << endl;
        }

        if(arguments.size())
        {
            QString command = arguments[0];

            if(!command.compare("-help" , Qt::CaseInsensitive))
            {
                cout << "-run [InputFilePath]" << endl;
                cout << "Example -run /users/curl/project.xml" << endl;
                cout << endl;
                cout << "-mrvmitem [MRVMItemType]" << endl;
                cout << "Produces the xml signature of MRVMItem xml" << endl;
                cout << "Example -mrvmitem RealMRVMItem" << endl;
                cout << endl;
                cout << "c" << endl;
                cout << "Close the program" << endl;

            }
            else if(!command.compare("-run" , Qt::CaseInsensitive))
            {
                if(arguments.size() == 2)
                {

                    QString filePath = arguments[1];

                    QFileInfo inputFile = QFileInfo(filePath);

                    if (QFile::exists(filePath))
                    {

                        cout << "Reading input file :" << filePath.toStdString() << endl;

                        MRVM mRVMModel(inputFile);

                        cout << "Finished reading input file :" << filePath.toStdString() << endl;

                        cout << "Running simulation..."<< endl;

                        mRVMModel.start();

                        cout << "Finished simulation..."<< endl;

                    }
                    else
                    {
                        cout << "A valid file path must be specified after -run" << endl;
                    }
                }
                else
                {
                    cout << "A valid file path must be specified after -run" << endl;
                }
            }
            else if(!command.compare("-mrvmitem"))
            {
                if(arguments.size() == 2)
                {
                    MRVMItem* item = nullptr;

                    if(!arguments[1].compare("RealMRVMItem" , Qt::CaseInsensitive))
                    {
                        item = new RealMRVMItem();
                    }
                    else if(!arguments[1].compare("RealArrayMRVMItem" , Qt::CaseInsensitive))
                    {
                        item = new RealArrayMRVMItem();
                    }
                    else if(!arguments[1].compare("CategoricalMRVMItem" , Qt::CaseInsensitive))
                    {
                        item = new CategoricalMRVMItem();
                        QMap<QString,int> test;
                        test["Category1"] = 1;
                        test["Category2"] = 2;
                        test["Category3"] = 3;

                        ((CategoricalMRVMItem*) item)->setCategories(test);
                    }
                    else if(!arguments[1].compare("RealRaster" , Qt::CaseInsensitive))
                    {
                        //item = new RealRaster();
                    }
                    else if(!arguments[1].compare("CategoricalRaster" , Qt::CaseInsensitive))
                    {
                        // item = new CategoricalRaster();
                        QMap<QString,int> test;
                        test["Category1"] = 1;
                        test["Category2"] = 2;
                        test["Category3"] = 3;

                        //((CategoricalMRVMItem*) item)->setCategories(test);
                    }

                    if(item)
                    {
                        cout << item->toString().toStdString() << endl;
                        delete item;
                    }
                    else
                    {
                        cout << "Specify a valid MRVMItem e.g., RealMRVMItem after -mrvmitem" << endl;
                    }


                }
                else
                {
                    cout << "Specify a valid MRVMItem e.g., RealMRVMItem after -mrvmitem" << endl;
                }

            }
            else if(!command.compare("c" , Qt::CaseInsensitive))
            {
                break;
            }
        }
    }

    return 0;
}



