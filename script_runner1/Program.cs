using System;
using System.Diagnostics;
using System.IO;
using System.Threading.Tasks;

namespace script_runner1
{
    class Program
    {
        static string window_masker(string path, string mpt_file, string x_min, string x_max, string y_min, string y_max)
        {
            string progToRun = "C:/Users/CJang/Desktop/Kyler_Speed_Circuit/python_notebooks/utils/window_masker.py";
            char[] splitter = { '\r' };

            Process proc = new Process();
            proc.StartInfo.FileName = "python.exe";
            proc.StartInfo.RedirectStandardOutput = true;
            proc.StartInfo.UseShellExecute = false;

            // call hello.py to concatenate passed parameters

            //"from tools import *; print(mpt_data(path=r'C:\Users\cjang\Desktop\Kyler_Speed_Circuit\data\\',data = ['DE_49_8_30.mpt']).guess_and_plot(mask = [1000018.6, 28]))"
            proc.StartInfo.Arguments = string.Concat(progToRun, " ", path, " ", mpt_file, " ", x_min, " ", x_max, " ", y_min, " ", y_max);
            //Console.WriteLine("Processing...");
            proc.Start();

            StreamReader sReader = proc.StandardOutput;
            string[] output = sReader.ReadToEnd().Split(splitter);

            //foreach (string s in output)
            //Console.WriteLine(s);
            string output1 = string.Join(", ", output);
            return output1;
        }
        static string mpt_plot(string path, string mpt_file)
        {
            string progToRun = "C:/Users/CJang/Desktop/Kyler_Speed_Circuit/python_notebooks/utils/mpt_plot.py";
            char[] splitter = { '\r' };

            Process proc = new Process();
            proc.StartInfo.FileName = "python.exe";
            proc.StartInfo.RedirectStandardOutput = true;
            proc.StartInfo.UseShellExecute = false;

            // call hello.py to concatenate passed parameters

            //"from tools import *; print(mpt_data(path=r'C:\Users\cjang\Desktop\Kyler_Speed_Circuit\data\\',data = ['DE_49_8_30.mpt']).guess_and_plot(mask = [1000018.6, 28]))"
            proc.StartInfo.Arguments = string.Concat(progToRun, " ", path, " ", mpt_file);
            //Console.WriteLine("Processing...");
            proc.Start();

            StreamReader sReader = proc.StandardOutput;
            string[] output = sReader.ReadToEnd().Split(splitter);

            //foreach (string s in output)
            //Console.WriteLine(s);
            string output1 = string.Join(", ", output);
            return output1;
        }
        static string masker(string path, string mpt_file, string mask_choice)
        {
            string progToRun = "C:/Users/CJang/Desktop/Kyler_Speed_Circuit/python_notebooks/utils/masker.py";
            char[] splitter = { '\r' };

            Process proc = new Process();
            proc.StartInfo.FileName = "python.exe";
            proc.StartInfo.RedirectStandardOutput = true;
            proc.StartInfo.UseShellExecute = false;

            // call hello.py to concatenate passed parameters

            //"from tools import *; print(mpt_data(path=r'C:\Users\cjang\Desktop\Kyler_Speed_Circuit\data\\',data = ['DE_49_8_30.mpt']).guess_and_plot(mask = [1000018.6, 28]))"
            proc.StartInfo.Arguments = string.Concat(progToRun, " ", path, " ", mpt_file , " ", mask_choice);
            Console.WriteLine("Processing...");
            proc.Start();

            StreamReader sReader = proc.StandardOutput;
            string[] output = sReader.ReadToEnd().Split(splitter);

            //foreach (string s in output)
            //Console.WriteLine(s);
            string output1 = string.Join(", ", output);
            return output1;
        }
        static string guesser(string path, string mpt_file)
        {
            string progToRun = "C:/Users/CJang/Desktop/Kyler_Speed_Circuit/python_notebooks/utils/guesser.py";
            char[] splitter = { '\r' };

            Process proc = new Process();
            proc.StartInfo.FileName = "python.exe";
            proc.StartInfo.RedirectStandardOutput = true;
            proc.StartInfo.UseShellExecute = false;

            // call hello.py to concatenate passed parameters

            //"from tools import *; print(mpt_data(path=r'C:\Users\cjang\Desktop\Kyler_Speed_Circuit\data\\',data = ['DE_49_8_30.mpt']).guess_and_plot(mask = [1000018.6, 28]))"
            proc.StartInfo.Arguments = string.Concat(progToRun, " ", path, " ", mpt_file);
            Console.WriteLine("Processing...");
            proc.Start();

            StreamReader sReader = proc.StandardOutput;
            string[] output = sReader.ReadToEnd().Split(splitter);

            //foreach (string s in output)
            //Console.WriteLine(s);
            string output1 = string.Join(", ", output);
            return output1;
        }
        static string path_listing(string path)
        {
            string progToRun = "C:/Users/CJang/Desktop/Kyler_Speed_Circuit/python_notebooks/utils/path_listing.py";
            char[] splitter = { '\r' };

            Process proc = new Process();
            proc.StartInfo.FileName = "python.exe";
            proc.StartInfo.RedirectStandardOutput = true;
            proc.StartInfo.UseShellExecute = false;

            // call hello.py to concatenate passed parameters

            //"from tools import *; print(mpt_data(path=r'C:\Users\cjang\Desktop\Kyler_Speed_Circuit\data\\',data = ['DE_49_8_30.mpt']).guess_and_plot(mask = [1000018.6, 28]))"
            proc.StartInfo.Arguments = string.Concat(progToRun, " ", path.ToString());
            proc.Start();

            StreamReader sReader = proc.StandardOutput;
            string[] output = sReader.ReadToEnd().Split(splitter);

            //foreach (string s in output)
            //Console.WriteLine(s);

            string output1 = string.Join(", ", output);
            if (output1.Length > 0)
            {
                return output1;
            }
            else
            {
                return ("Empty folder or a bad directory");
            }
        }
        static string mpt_dataframe(string path, string mpt_file)
        {
            Process proc = new Process();
            proc.StartInfo.FileName = "python.exe";
            proc.StartInfo.RedirectStandardOutput = true;
            proc.StartInfo.UseShellExecute = false;

            string progToRun = "C:/Users/CJang/Desktop/Kyler_Speed_Circuit/python_notebooks/utils/mpt_dataframe.py";
            char[] splitter = { '\r' };

            // call hello.py to concatenate passed parameters

            //"from tools import *; print(mpt_data(path=r'C:\Users\cjang\Desktop\Kyler_Speed_Circuit\data\\',data = ['DE_49_8_30.mpt']).guess_and_plot(mask = [1000018.6, 28]))"
            proc.StartInfo.Arguments = string.Concat(progToRun, " ", path, " ", mpt_file);
            proc.Start();

            StreamReader sReader = proc.StandardOutput;
            string[] output = sReader.ReadToEnd().Split(splitter);

            string output1 = string.Join(", ", output);
            return output1;
        }
        static void Main(string[] args)
        {
            string path;
            string data;
            string z;
            string x_min;
            string x_max;
            string y_min;
            string y_max;



            Console.WriteLine("Pick a path: ");
            try
            {
                path = (Console.ReadLine());
            }
            catch (Exception e)
            {
                Console.WriteLine("Not a String; Pick again");
                path = (Console.ReadLine());
            }

            Console.WriteLine(path_listing(path));


            if (path_listing(path) != "Empty folder or a bad directory")
            {
                Console.WriteLine("Pick a datafile from the file above: ");
                try
                {
                    data = (Console.ReadLine());
                }
                catch (Exception e)
                {
                    Console.WriteLine("Not a String; Pick again");
                    data = (Console.ReadLine());
                }

                Console.WriteLine(guesser(path, data));


                Console.WriteLine("Now pick a mask: ");
                Console.WriteLine("1: Fast Mask");
                Console.WriteLine("2: Binning Mask");
                Console.WriteLine("3: Combo Mask");

                try
                {
                    z = (Console.ReadLine());
                }
                catch (Exception e)
                {
                    Console.WriteLine("Not applicable");
                    z = (Console.ReadLine());
                }

                Console.WriteLine(masker(path, data, z));


            }
            else
            {
                Console.WriteLine("Bad path");
            }
            



            

            //Console.WriteLine("Pick X Limits: ");
            //try
            //{
            //    x_min = (Console.ReadLine());
            //    x_max = (Console.ReadLine());
            //}
            //catch (Exception e)
            //{
            //   Console.WriteLine("Not a String; Pick again");
            //    x_min = (Console.ReadLine());
            //    x_max = (Console.ReadLine());
            //}
            
            //Console.WriteLine("Pick Y Limits: ");
            //try
            //{
            //    y_min = (Console.ReadLine());
            //    y_max = (Console.ReadLine());
            //}
            //catch (Exception e)
            //{
            //    Console.WriteLine("Not a String; Pick again");
            //    y_min = (Console.ReadLine());
            //   y_max = (Console.ReadLine());
            //}

            //Console.WriteLine(window_masker(path, data, x_min, x_max, y_min, y_max));


        }
    }
}