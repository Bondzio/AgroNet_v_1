﻿//------------------------------------------------------------------------------
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
    using System.Data.Entity;
    using System.Data.Entity.Infrastructure;
    
    public partial class PLMUsers : DbContext
    {
        public PLMUsers()
            : base("name=PLMUsers")
        {
        }
    
        protected override void OnModelCreating(DbModelBuilder modelBuilder)
        {
            throw new UnintentionalCodeFirstException();
        }
    
        public DbSet<Users> Users { get; set; }
        public DbSet<ActivitySessions> ActivitySessions { get; set; }
        public DbSet<Applications> Applications { get; set; }
        public DbSet<ApplicationUsers> ApplicationUsers { get; set; }
        public DbSet<WebPages> WebPages { get; set; }
        public DbSet<Roles> Roles { get; set; }
        public DbSet<ActivityLogs> ActivityLogs { get; set; }
        public DbSet<OperationRoles> OperationRoles { get; set; }
        public DbSet<Operations> Operations { get; set; }
        public DbSet<Tables> Tables { get; set; }
        public DbSet<UserCountries> UserCountries { get; set; }
        public DbSet<CountriesUser> CountriesUser { get; set; }
    }
}
