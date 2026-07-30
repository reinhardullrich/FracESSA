using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Diagnostics;
using System.IO;

namespace EFR
{
    class EssFinderQueue
    {
        readonly object _locker = new object();

        private Thread[] _workers;
        private Queue<EssFinder> _itemQ = new Queue<EssFinder>();
        StreamWriter _tempm = new StreamWriter(MainClass.tempm);
        StreamWriter _tempv = new StreamWriter(MainClass.tempv);

        public EssFinderQueue()
        {
            _workers = new Thread[MainClass.NumberOfWorkerThreads];

            StringBuilder hm = new StringBuilder("");
            hm.Append("'matrixid',");
            hm.Append("'matrixname',");
            hm.Append("'matrix',");
            hm.Append("'dimension',");
            hm.Append("'anzahless'");
            _tempm.WriteLine(hm.ToString());

            StringBuilder hv = new StringBuilder("");
            hv.Append("'matrixid',");
            hv.Append("'vektor',");
            hv.Append("'supportsize',");
            hv.Append("'extsupportsize',");
            hv.Append("'support',");
            hv.Append("'extsupport',");
            hv.Append("'isess',");
            hv.Append("'grundess',");
            hv.Append("'payoff'");
            _tempv.WriteLine(hv.ToString());

            // Create and start a separate thread for each worker
            for (int i = 0; i < MainClass.NumberOfWorkerThreads; i++)
            {
                (_workers[i] = new Thread(Consume)).Start();
                _workers[i].Name = "Worker " + i;
            }
        }

        
        public void EnqueueItem(EssFinder item)
        {
            lock (_locker)
            {
                while (_itemQ.Count > MainClass.MaxNumberPCQueue)
                    Monitor.Wait(_locker, MainClass.PCEnqueueWaitingTime);
                _itemQ.Enqueue(item);
                Monitor.Pulse(_locker); // We must pulse because we're changing a blocking condition.
            }
        }

        private void Consume()
        {
            while (true) // Keep consuming until told otherwise.
            {
                EssFinder item;
                lock (_locker)
                {
                    while (_itemQ.Count == 0) 
                        Monitor.Wait (_locker);
                    item = _itemQ.Dequeue();                    
                }
                if (item == null) 
                    return; // This signals our exit.

                // Execute item.
                item.FindESS();
                //Console.WriteLine(item.MatrixId);

                StringBuilder strm = new StringBuilder("");
                StringBuilder strv = new StringBuilder("");

                strm.Append("'" + item.MatrixId + "',");
                strm.Append("'" + item.Name + "',");
                strm.Append("'',"); // strm.Append("'" + item.SpielMatrix.ToRString() + "',"); // strm.Append("'',"); // 
                strm.Append("'" + item.Dimension + "',");
                strm.Append("'" + item.AnzahlEss + "'");
                
                foreach (VectorInformation vector in item.VectorSet)
                {  
                    strv.Append("'" + item.MatrixId + "',");
                    strv.Append("'',"); // strv.Append("'" + vector.Vektor + "',"); //strv.Append("'',"); 
                    strv.Append("'" + vector.SupportSize + "',");
                    strv.Append("'" + vector.ExtendedSupportSize + "',");
                    strv.Append("'" + vector.Support + "',");
                    strv.Append("'" + (vector.Support != vector.ExtendedSupport ?  vector.ExtendedSupport : "") + "',");
                    strv.Append("'" + (vector.IsEss == true ? 1 : 0) + "',");
                    strv.Append("'" + vector.GrundEss + "',");
                    strv.Append("'" + vector.Payoff.ToString().Replace(",", ".") + "'");
                    strv.Append(Environment.NewLine);                   
                }

                lock (_locker)
                {
                    _tempm.WriteLine(strm.ToString());
                    _tempv.Write(strv.ToString());
                    //_tempm.Flush();
                    //_tempv.Flush();

                }
            }
        }

        public void Shutdown()
        {           
            // Enqueue one null item per worker to make each exit.
            foreach (Thread worker in _workers)
                EnqueueItem(null);
            // Wait for workers to finish
            foreach (Thread worker in _workers)
                worker.Join();
            _tempm.Flush();
            _tempm.Close();
            _tempm.Dispose();
            _tempv.Flush();
            _tempv.Close();
            _tempv.Dispose();
        }

    }
}
