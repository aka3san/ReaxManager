using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ReaxManager
{
    public class Molecule
    {
        private int id;
        private string smiles;
        public Molecule(int id, string smiles)
        {
            this.id = id;
            this.smiles = smiles;
        }
    }
}
