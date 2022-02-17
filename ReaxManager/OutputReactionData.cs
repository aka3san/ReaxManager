using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
namespace ReaxManager
{
    public class OutputReactionData
    {
        private int totalTime;
        private List<int> totalMolNum;

        public OutputReactionData(int totalTime)
        {
            this.totalTime = totalTime;
        }

        public void OutputToTextFile()
        {
            try
            {
                File.WriteAllText(@"ReactionData.txt", "Reaction Data" + Environment.NewLine);
            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of Molecules" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalMolNum[i]}" + Environment.NewLine);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of D4OH" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalD4OHNum[i]}" + Environment.NewLine);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of H2O" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalH2ONum[i]}" + Environment.NewLine);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of O2" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalO2Num[i]}" + Environment.NewLine);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of Reaction" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalReactionNum[i]}" + Environment.NewLine);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of Species" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalSpecies[i]}" + Environment.NewLine);
            }


            List<int> polymolCountPerTime = new List<int>();
            File.AppendAllText(@"ReactionData.txt", "List of Species" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                //重合した分子の数
                int polymolCount = 0;

                File.AppendAllText(@"ReactionData.txt", $"TimeStep{i}" + Environment.NewLine);
                foreach (KeyValuePair<string, int> smiles in totalSpeciesDict[i])
                {
                    string molFormula = "";
                    int cCount, hCount, fCount, oCount;
                    cCount = 0; hCount = 0; fCount = 0; oCount = 0;
                    foreach (char atomName in totalSmilesPlusH[i][smiles.Key])
                    {
                        switch (atomName)
                        {
                            case 'C':
                                cCount++;
                                break;
                            case 'H':
                                hCount++;
                                break;
                            case 'F':
                                fCount++;
                                break;
                            case 'O':
                                oCount++;
                                break;
                        }
                    }
                    if (cCount > 37)
                    {
                        polymolCount++;
                    }
                    molFormula += (cCount != 0) ? ("C" + cCount.ToString()) : "";
                    molFormula += (hCount != 0) ? ("H" + hCount.ToString()) : "";
                    molFormula += (fCount != 0) ? ("F" + fCount.ToString()) : "";
                    molFormula += (oCount != 0) ? ("O" + oCount.ToString()) : "";
                    if (i == 0 || i == totalTime - 2)
                    {
                        File.AppendAllText(@"ReactionData.txt", $"●{smiles.Key} \"({molFormula})\": {smiles.Value}" + Environment.NewLine);
                    }
                }
                polymolCountPerTime.Add(polymolCount);
            }

            File.AppendAllText(@"ReactionData.txt", $"●重合した分子の数" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{polymolCountPerTime[i]}" + Environment.NewLine);
            }

            File.AppendAllText(@"ReactionData.txt", "ReactionPerTime" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"TimeStep{i}" + Environment.NewLine);
                for (int i2 = 0; i2 < totalReaction[i].Count; i2++)
                {
                    File.AppendAllText(@"ReactionData.txt", $"------反応{i2 + 1}------" + Environment.NewLine);
                    for (int j = 0; j < totalReaction[i][i2].Count; j++)
                    {
                        for (int k = 0; k < totalReaction[i][i2][j].Count; k++)
                        {
                            File.AppendAllText(@"ReactionData.txt", $"{totalReaction[i][i2][j][k]}" + Environment.NewLine);
                            if (k == totalReaction[i][i2][j].Count - 1)
                            {
                                break;
                            }
                            File.AppendAllText(@"ReactionData.txt", "+" + Environment.NewLine);
                        }
                        if (j == 0)
                        {
                            File.AppendAllText(@"ReactionData.txt", "--↓↓--" + Environment.NewLine);
                        }
                    }
                    File.AppendAllText(@"ReactionData.txt", "------------" + Environment.NewLine);
                }
            }
        }
    }
}
