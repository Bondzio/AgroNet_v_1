//------------------------------------------------------------------------------
// <auto-generated>
//    Este código se generó a partir de una plantilla.
//
//    Los cambios manuales en este archivo pueden causar un comportamiento inesperado de la aplicación.
//    Los cambios manuales en este archivo se sobrescribirán si se regenera el código.
// </auto-generated>
//------------------------------------------------------------------------------

namespace AgroNet.Models
{
    using System;
    using System.Collections.Generic;
    
    public partial class DivisionImagesSizes
    {
        public int DivisionImageId { get; set; }
        public byte ImageSizeId { get; set; }
    
        public virtual DivisionImages DivisionImages { get; set; }
        public virtual ImageSizes ImageSizes { get; set; }
    }
}
